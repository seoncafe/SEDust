module mc_engine
   ! Thread-safe MC engine.  Takes an mc_grain_t (filled by one of the
   ! grain_setup_* routines) and an rng_t and produces a P(T) histogram
   ! plus aggregate (t, T) records.  Each thread should own its own
   ! mc_grain_t and rng_t.  The only shared module state is the optional
   ! trajectory buffer (single-thread use) and the step-cap diagnostic
   ! counters, which are updated inside an OpenMP critical region.

   use constants,      only: wp, pi, c_cgs => c, h_cgs => h_planck, c_si, k_boltz
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use mc_grain_type,  only: mc_grain_t
   use mc_rng,         only: rng_t, rng_uniform, rng_exp
   use mc_heatcap,     only: U_of_T, C_of_T
   use radfield,       only: J_Mathis, bbody
   use q_graphite_d16_sphere_mod, only: q_graphite_d16_sphere_abs
   implicit none
   private
   public :: hc_erg_um, sample_photon_lam
   public :: grain_setup_graphite, grain_setup_from_cabs
   public :: grain_get_C, grain_get_QPL, grain_T_from_U
   public :: grain_radiative_equilibrium, mc_equilibrate
   public :: mc_run_engine
   public :: mc_run_engine_2pass, mc_run_engine_buffered
   public :: mc_traj_enable, mc_traj_disable, mc_traj_get
   public :: MC_OK, MC_STEP_BUDGET, MC_INVALID_STEP, MC_INVALID_WAIT
   public :: MC_HISTOGRAM_RANGE, MC_SAMPLE_CAPACITY

   integer, parameter :: MC_OK = 0, MC_STEP_BUDGET = 1
   integer, parameter :: MC_INVALID_STEP = 2, MC_INVALID_WAIT = 3
   integer, parameter :: MC_HISTOGRAM_RANGE = 4, MC_SAMPLE_CAPACITY = 5
   integer, parameter :: DEFAULT_STEP_BUDGET = 1000000

   real(wp), parameter :: hc_erg_um   = h_cgs * c_cgs * 1.0e4_wp
   real(wp), parameter :: sigma_sb_cgs= 2.0_wp*pi**5*k_boltz**4/(15.0_wp*h_cgs**3*c_cgs**2)

   ! ------------------------------------------------------------------
   ! Trajectory dump: when enabled via mc_traj_enable(stride, cap), the
   ! cool_segment routine pushes (t_global, T) at every adaptive sub-step
   ! into a module-level buffer, subject to a minimum time stride.  The
   ! event boundaries (T_pre, T_post) are pushed separately by
   ! mc_run_engine.  After the run, main_mc calls mc_traj_get to retrieve
   ! the buffer and writes it to _Tt.dat.  Stride sampling keeps the
   ! output bounded for runs with many events; for visualization use only.
   ! ------------------------------------------------------------------
   logical,            save :: traj_enabled = .false.
   logical,            save :: traj_overflow = .false.
   integer,            save :: traj_n       = 0
   integer,            save :: traj_cap     = 0
   real(wp),           save :: traj_stride  = 0.0_wp
   real(wp),           save :: traj_t_last  = -1.0e30_wp
   ! Set by mc_run_engine before each cool_segment call so that
   ! cool_segment can map its internal tau_seg to a global wall time.
   real(wp),           save :: traj_seg_start_time = 0.0_wp
   real(wp), allocatable, save :: traj_time_buf(:), traj_temp_buf(:)

   ! ------------------------------------------------------------------
   ! Adaptive-engine diagnostics count completed cooling segments and
   ! segments needing more than one chunk of integration steps.
   ! ------------------------------------------------------------------
   integer(kind=8), save :: stepcap_nseg = 0_8
   integer(kind=8), save :: stepcap_nhit = 0_8

contains

   subroutine mc_traj_enable(stride_in, cap_in)
      ! Allocate and turn on the trajectory buffer.
      real(wp), intent(in) :: stride_in
      integer,  intent(in) :: cap_in
      if (allocated(traj_time_buf)) deallocate(traj_time_buf, traj_temp_buf)
      allocate(traj_time_buf(cap_in), traj_temp_buf(cap_in))
      traj_cap     = cap_in
      traj_stride  = stride_in
      traj_n       = 0
      traj_t_last  = -1.0e30_wp
      traj_enabled = .true.
      traj_overflow = .false.
   end subroutine mc_traj_enable


   subroutine mc_traj_disable()
      if (allocated(traj_time_buf)) deallocate(traj_time_buf, traj_temp_buf)
      traj_enabled = .false.
      traj_overflow = .false.
      traj_n = 0
      traj_cap = 0
   end subroutine mc_traj_disable


   subroutine mc_traj_get(n_out, time_out, temp_out)
      ! Copy the buffer contents to caller arrays.  Each output array must
      ! be at least traj_n long; pass size(.)=traj_cap to be safe.
      integer,  intent(out) :: n_out
      real(wp), intent(out) :: time_out(:), temp_out(:)
      n_out = traj_n
      if (n_out > 0) then
         time_out(1:n_out) = traj_time_buf(1:n_out)
         temp_out(1:n_out) = traj_temp_buf(1:n_out)
      end if
   end subroutine mc_traj_get


   subroutine traj_push(t_global, T_now, force)
      ! Append (t_global, T_now) to the buffer if (a) buffer is enabled,
      ! (b) buffer has room, and (c) either force=.true. (e.g. event
      ! boundary) or the elapsed time since the last record exceeds the
      ! stride.
      real(wp), intent(in) :: t_global, T_now
      logical,  intent(in) :: force
      if (.not. traj_enabled) return
      if (.not. force) then
         if (t_global < traj_t_last + traj_stride) return
      end if
      if (traj_n >= traj_cap) then
         traj_overflow = .true.
         return
      end if
      traj_n = traj_n + 1
      traj_time_buf(traj_n) = t_global
      traj_temp_buf(traj_n) = T_now
      traj_t_last    = t_global
   end subroutine traj_push


   subroutine stepcap_tally(tag, a_um, nseg_call, ncap_call)
      ! Add this adaptive-engine call's cooling-segment and chunk-continuation
      ! counts to the process-wide totals and print a one-line diagnostic
      ! summary.  Called once per engine invocation (not per step), inside a
      ! named OpenMP critical region so the size loop in mc_sed stays
      ! thread-safe.  The cumulative fields let a reader recover the run total
      ! from the last line (e.g. tail -1).
      character(len=*), intent(in) :: tag
      real(wp),         intent(in) :: a_um
      integer(kind=8),  intent(in) :: nseg_call, ncap_call
      !$omp critical (stepcap_report)
      stepcap_nseg = stepcap_nseg + nseg_call
      stepcap_nhit = stepcap_nhit + ncap_call
      write(*,'(a,a,a,es10.2,a,i0,a,i0,a,i0,a,i0)') &
         ' [stepcap ', trim(tag), ' a=', a_um, &
         ' um] seg=', nseg_call, ' chunk_continuations=', ncap_call, &
         ' | cum seg=', stepcap_nseg, ' chunk_continuations=', stepcap_nhit
      !$omp end critical (stepcap_report)
   end subroutine stepcap_tally

   ! =========================================================================
   ! Setup variants
   ! =========================================================================

   subroutine grain_setup_graphite(g, a_um, U_isrf, comp, lam_c_in, NLAM_in)
      ! Set up a grain using the D16 graphite sphere Q from
      ! q_graphite_d16_sphere_mod.  This is the standalone-validation path
      ! used for PIIM Figures 24.5 / 24.6.
      type(mc_grain_t), intent(out)          :: g
      real(wp),         intent(in)           :: a_um, U_isrf
      character(len=*), intent(in)           :: comp
      real(wp),         intent(in), optional :: lam_c_in
      integer,          intent(in), optional :: NLAM_in
      real(wp), allocatable :: Cabs_tmp(:), lam_tmp(:)
      integer  :: NL, i
      real(wp) :: lam_min, lam_max, dlnlam

      NL = 2000
      if (present(NLAM_in)) NL = NLAM_in
      if (NL < 2) error stop 'grain_setup_graphite: at least two wavelengths required'

      allocate(lam_tmp(NL), Cabs_tmp(NL))
      lam_min = 1.0e-3_wp
      lam_max = 1.0e4_wp
      dlnlam  = log(lam_max/lam_min) / real(NL-1, wp)
      do i = 1, NL
         lam_tmp(i) = lam_min * exp((i-1)*dlnlam)
         call q_graphite_d16_sphere_abs(a_um, lam_tmp(i), Cabs_tmp(i))
         if (Cabs_tmp(i) < 0.0_wp) Cabs_tmp(i) = 0.0_wp
      end do

      call grain_setup_from_cabs(g, a_um, U_isrf, comp, lam_tmp, Cabs_tmp, &
                                 lam_c_in=lam_c_in)
      deallocate(lam_tmp, Cabs_tmp)
   end subroutine grain_setup_graphite


   subroutine grain_setup_from_cabs(g, a_um, U_isrf, comp, lam_in, Q_in, lam_c_in, J_lam_in)
      ! General setup: caller supplies the wavelength grid and Q_abs(lam).
      ! Used by the SED builder for arbitrary grain compositions
      ! (e.g., astrodust grains, PAH grains with xi-blend).
      type(mc_grain_t), intent(out)          :: g
      real(wp),         intent(in)           :: a_um, U_isrf
      character(len=*), intent(in)           :: comp
      real(wp),         intent(in)           :: lam_in(:), Q_in(:)
      real(wp),         intent(in), optional :: lam_c_in
      real(wp),         intent(in), optional :: J_lam_in(:) ! W m^-3 sr^-1

      integer  :: NL, NCDF, NT, i, j
      real(wp), allocatable :: J_SI(:), integrand_cont(:), integrand_event(:)
      real(wp), allocatable :: cum_event(:), BB(:), lam_node(:)
      real(wp) :: w, Fmax, integ_BQ, integ_B, Tt, frac_c, ic_c
      real(wp) :: photon_min, q_edge, field_edge(1), edge(1), frac_edge

      NL = size(lam_in)
      if (NL < 2) error stop 'grain_setup_from_cabs: at least two wavelengths required'
      if (.not. all(ieee_is_finite(lam_in)) .or. any(lam_in <= 0.0_wp)) &
         error stop 'grain_setup_from_cabs: wavelengths must be finite and positive'
      if (any(lam_in(2:) <= lam_in(:NL-1))) &
         error stop 'grain_setup_from_cabs: wavelengths must be increasing'
      if (.not. ieee_is_finite(a_um) .or. a_um <= 0.0_wp) &
         error stop 'grain_setup_from_cabs: radius must be finite and positive'
      if (size(Q_in) /= NL) then
         write(*,'(a)') 'grain_setup_from_cabs: size(lam_in) /= size(Q_in)'
         stop 1
      end if
      if (present(J_lam_in)) then
         if (size(J_lam_in) /= NL) error stop 'grain_setup_from_cabs: J_lam_in size mismatch'
         if (.not. all(ieee_is_finite(J_lam_in))) &
            error stop 'grain_setup_from_cabs: J_lam_in must be finite'
         if (any(J_lam_in < 0.0_wp)) &
            error stop 'grain_setup_from_cabs: J_lam_in must be nonnegative'
      end if

      if (size(Q_in) /= NL) error stop 'grain_setup_from_cabs: opacity size mismatch'
      if (.not. all(ieee_is_finite(Q_in)) .or. any(Q_in < 0.0_wp)) &
         error stop 'grain_setup_from_cabs: opacity must be finite and nonnegative'
      g%a_um = a_um
      g%comp = comp
      if (present(lam_c_in)) g%lam_c = lam_c_in

      if (.not. ieee_is_finite(g%lam_c) .or. g%lam_c <= 0.0_wp) &
         error stop 'grain_setup_from_cabs: cutoff must be finite and positive'
      g%NLAM = NL
      allocate(g%lam_grid(NL), g%Q_grid(NL), g%u_lam_grid(NL))
      g%lam_grid = lam_in
      g%Q_grid   = Q_in

      ! Radiation field u_lam from the supplied field or legacy Mathis.
      !   u_lam_cgs[erg/cm3/um] = (4 pi / c_SI) * J_SI[W/m2/m/sr] * 1e-5
      allocate(J_SI(NL))
      if (present(J_lam_in)) then
         J_SI = J_lam_in
      else
         call J_Mathis(U_isrf, g%lam_grid, J_SI)
      end if
      g%u_lam_grid = (4.0_wp * pi / c_si) * J_SI * 1.0e-5_wp
      deallocate(J_SI)

      ! Continuous + stochastic split
      allocate(integrand_cont(NL), integrand_event(NL))
      do i = 1, NL
         integrand_cont(i)  = g%Q_grid(i) * g%u_lam_grid(i)
         integrand_event(i) = pi * (a_um*1.0e-4_wp)**2 * g%Q_grid(i) &
                              * c_cgs * g%u_lam_grid(i) * g%lam_grid(i) / hc_erg_um
      end do
      g%H_cont = 0.0_wp
      do i = 1, NL-1
         w = 0.5_wp * (g%lam_grid(i+1)-g%lam_grid(i))
         if (g%lam_grid(i+1) <= g%lam_c) then
            cycle
         else if (g%lam_grid(i) >= g%lam_c) then
            g%H_cont = g%H_cont+w*(integrand_cont(i)+integrand_cont(i+1))
         else
            frac_c = (g%lam_c-g%lam_grid(i))/(g%lam_grid(i+1)-g%lam_grid(i))
            ic_c = integrand_cont(i)+frac_c*(integrand_cont(i+1)-integrand_cont(i))
            g%H_cont = g%H_cont+0.5_wp*(g%lam_grid(i+1)-g%lam_c)*(ic_c+integrand_cont(i+1))
         end if
      end do
      g%H_cont = (c_cgs / 4.0_wp) * g%H_cont

      ! Store the CDF on the physical wavelength knots, not a uniformly
      ! sampled inverse. Interpolating a sparse inverse smears its jumps
      ! across intervals with zero absorption probability. The density is
      ! nonnegative and linear within each original wavelength interval.
      NCDF = NL
      g%NCDF = NCDF
      allocate(g%cdf_F(NCDF), g%cdf_lam(NCDF), g%cdf_density(NCDF))
      allocate(cum_event(NL), lam_node(NL))
      photon_min = g%lam_grid(1)
      if (.not. present(J_lam_in)) then
         ! The legacy Mathis stellar field is zero below the Lyman edge.
         ! Clip only photon integration, retaining the full cooling grid.
         photon_min = max(photon_min, 0.0912_wp)
      end if
      do i = 1, NL
         lam_node(i) = min(max(g%lam_grid(i), photon_min), g%lam_c)
         if (lam_node(i) == g%lam_grid(i)) then
            g%cdf_density(i) = integrand_event(i)
         else
            j = 1
            do while (j < NL-1 .and. g%lam_grid(j+1) < lam_node(i))
               j = j + 1
            end do
            frac_edge = (lam_node(i)-g%lam_grid(j))/(g%lam_grid(j+1)-g%lam_grid(j))
            g%cdf_density(i) = integrand_event(j) + frac_edge*(integrand_event(j+1)-integrand_event(j))
            if (.not. present(J_lam_in)) then
               ! Evaluate the discontinuous field at its supported edge;
               ! interpolation through a zero UV node would invent a ramp.
               q_edge = g%Q_grid(j) + frac_edge*(g%Q_grid(j+1)-g%Q_grid(j))
               edge(1) = lam_node(i)
               call J_Mathis(U_isrf, edge, field_edge)
               g%cdf_density(i) = pi*(a_um*1.0e-4_wp)**2 * q_edge &
                  * (4.0_wp*pi/c_si)*field_edge(1)*1.0e-5_wp*c_cgs*lam_node(i)/hc_erg_um
            end if
         end if
      end do
      cum_event(1) = 0.0_wp
      do i = 2, NL
         cum_event(i) = cum_event(i-1) + 0.5_wp*(lam_node(i)-lam_node(i-1)) &
            * (g%cdf_density(i-1)+g%cdf_density(i))
      end do
      Fmax = cum_event(NL)
      if (.not. ieee_is_finite(Fmax) .or. Fmax <= 0.0_wp) &
         error stop 'grain_setup_from_cabs: no supported stochastic absorption'
      g%rate_event = Fmax
      g%cdf_F = cum_event/Fmax
      g%cdf_lam = lam_node
      g%cdf_density = g%cdf_density/Fmax
      deallocate(cum_event, lam_node, integrand_cont, integrand_event)

      ! Thermo tables
      NT = 2000
      g%NT = NT
      allocate(g%T_grid(NT), g%U_grid(NT), g%QPL_grid(NT))
      do i = 1, NT
         g%T_grid(i) = 1.0_wp * (1.0e4_wp/1.0_wp)**(real(i-1,wp)/real(NT-1,wp))
      end do
      do i = 1, NT
         g%U_grid(i) = U_of_T(g%T_grid(i), a_um, comp)
      end do
      allocate(BB(NL))
      do i = 1, NT
         Tt = g%T_grid(i)
         do j = 1, NL
            BB(j) = bbody(Tt, g%lam_grid(j))
         end do
         integ_BQ = 0.0_wp
         integ_B  = 0.0_wp
         do j = 1, NL-1
            w = 0.5_wp * (g%lam_grid(j+1) - g%lam_grid(j))
            integ_BQ = integ_BQ + w * (g%Q_grid(j)*BB(j) + g%Q_grid(j+1)*BB(j+1))
            integ_B  = integ_B  + w * (BB(j) + BB(j+1))
         end do
         if (integ_B > 0.0_wp) then
            g%QPL_grid(i) = integ_BQ / integ_B
         else
            g%QPL_grid(i) = 0.0_wp
         end if
      end do
      deallocate(BB)
   end subroutine grain_setup_from_cabs


   ! =========================================================================
   ! Accessors
   ! =========================================================================

   function grain_get_C(g, Tval) result(C)
      type(mc_grain_t), intent(in) :: g
      real(wp),         intent(in) :: Tval
      real(wp) :: C
      C = C_of_T(Tval, g%a_um, g%comp)
   end function grain_get_C


   function grain_get_QPL(g, Tval) result(Q)
      ! Planck-averaged <Q>_T from the pre-tabulated grid by log-T interpolation.
      ! Outside the table, hold its endpoint value as a numerical continuation.
      ! This does not resolve the physical Planck mean below 1 K. Its use for
      ! emission requires negligible power from that unresolved temperature tail.
      type(mc_grain_t), intent(in) :: g
      real(wp),         intent(in) :: Tval
      real(wp) :: Q, lT, lT1, lT2, frac
      integer  :: lo, hi, mid
      if (Tval <= g%T_grid(1)) then;  Q = g%QPL_grid(1);  return; end if
      if (Tval >= g%T_grid(g%NT)) then; Q = g%QPL_grid(g%NT); return; end if
      lT = log(Tval)
      lo = 1; hi = g%NT
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (g%T_grid(mid) <= Tval) then;  lo = mid;  else;  hi = mid;  end if
      end do
      lT1 = log(g%T_grid(lo));  lT2 = log(g%T_grid(hi))
      frac = (lT - lT1) / (lT2 - lT1)
      Q = g%QPL_grid(lo) + frac * (g%QPL_grid(hi) - g%QPL_grid(lo))
   end function grain_get_QPL


   function grain_T_from_U(g, U_target, status_out) result(T)
      ! Invert the analytic enthalpy, using the table only to locate a bracket.
      ! The 0.5 K numerical floor is the lowest representable temperature.
      ! Above the table, expand the bracket rather than capping the energy.
      type(mc_grain_t), intent(in) :: g
      real(wp), intent(in) :: U_target
      integer, intent(out), optional :: status_out
      real(wp) :: T, T_lo, T_hi, U_lo, U_hi, U_now, C_now, T_new, tol
      integer :: lo, hi, mid, iter, status
      integer, parameter :: MAX_INVERT_STEPS = 64
      status = MC_OK
      T = 0.0_wp
      T_lo = 0.0_wp
      T_hi = 0.0_wp
      U_lo = 0.0_wp
      U_hi = 0.0_wp
      if (.not. ieee_is_finite(U_target) .or. U_target <= 0.0_wp) then
         status = MC_INVALID_STEP
      else
         if (U_target <= g%U_grid(1)) then
            T_lo = 0.5_wp
            U_lo = U_of_T(T_lo, g%a_um, g%comp)
            if (.not. ieee_is_finite(U_lo) .or. U_lo <= 0.0_wp) then
               status = MC_INVALID_STEP
            else if (U_target < U_lo) then
               status = MC_INVALID_STEP
            else if (U_target <= U_lo) then
               T = T_lo
            else
               T_hi = g%T_grid(1)
               U_hi = g%U_grid(1)
            end if
         else if (U_target <= g%U_grid(g%NT)) then
            lo = 1
            hi = g%NT
            do while (hi - lo > 1)
               mid = (lo + hi) / 2
               if (g%U_grid(mid) <= U_target) then
                  lo = mid
               else
                  hi = mid
               end if
            end do
            T_lo = g%T_grid(lo)
            T_hi = g%T_grid(hi)
            U_lo = g%U_grid(lo)
            U_hi = g%U_grid(hi)
         else
            T_lo = g%T_grid(g%NT)
            U_lo = g%U_grid(g%NT)
            T_hi = T_lo
            U_hi = U_lo
            if (U_target > U_hi) then
               do iter = 1, MAX_INVERT_STEPS
                  if (T_hi > huge(T_hi) / 2.0_wp) exit
                  T_hi = 2.0_wp * T_hi
                  U_hi = U_of_T(T_hi, g%a_um, g%comp)
                  if (.not. ieee_is_finite(U_hi)) exit
                  if (U_hi >= U_target) exit
                  T_lo = T_hi
                  U_lo = U_hi
               end do
            end if
         end if

         if (status == MC_OK .and. T <= 0.0_wp) then
            if (.not. ieee_is_finite(U_hi) .or. U_hi < U_target .or. &
                U_hi <= U_lo .or. T_hi <= T_lo) then
               status = MC_INVALID_STEP
            else
               tol = 64.0_wp * epsilon(U_target) * max(abs(U_target), abs(U_hi))
               T = T_lo + (T_hi - T_lo) * (U_target - U_lo) / (U_hi - U_lo)
               T = max(T_lo, min(T_hi, T))
               do iter = 1, MAX_INVERT_STEPS
                  U_now = U_of_T(T, g%a_um, g%comp)
                  if (.not. ieee_is_finite(U_now)) exit
                  if (abs(U_now - U_target) <= tol) exit
                  if (U_now < U_target) then
                     T_lo = T
                     U_lo = U_now
                  else
                     T_hi = T
                     U_hi = U_now
                  end if
                  C_now = C_of_T(T, g%a_um, g%comp)
                  if (ieee_is_finite(C_now) .and. C_now > 0.0_wp) then
                     T_new = T - (U_now - U_target) / C_now
                  else
                     T_new = 0.5_wp * (T_lo + T_hi)
                  end if
                  if (.not. ieee_is_finite(T_new)) T_new = 0.5_wp * (T_lo + T_hi)
                  if (T_new <= T_lo .or. T_new >= T_hi .or. (T_new <= T .and. T_new >= T)) &
                     T_new = 0.5_wp * (T_lo + T_hi)
                  if ((T_new <= T .and. T_new >= T)) exit
                  T = T_new
               end do
               U_now = U_of_T(T, g%a_um, g%comp)
               if (.not. ieee_is_finite(U_now)) then
                  status = MC_INVALID_STEP
               else if (abs(U_now - U_target) > tol) then
                  ! When adjacent representable temperatures bracket the
                  ! target, analytic U(T) can jump across it numerically.
                  ! Return the endpoint with the smaller enthalpy residual.
                  if (T_hi <= nearest(T_lo, 1.0_wp)) then
                     U_lo = U_of_T(T_lo, g%a_um, g%comp)
                     U_hi = U_of_T(T_hi, g%a_um, g%comp)
                     if (ieee_is_finite(U_lo) .and. ieee_is_finite(U_hi) .and. &
                         U_target >= min(U_lo,U_hi) .and. U_target <= max(U_lo,U_hi) .and. &
                         min(abs(U_lo-U_target),abs(U_hi-U_target)) <= &
                         1.0e-12_wp * max(abs(U_target),abs(U_lo),abs(U_hi))) then
                        if (abs(U_lo - U_target) <= abs(U_hi - U_target)) then
                           T = T_lo
                        else
                           T = T_hi
                        end if
                     else
                        status = MC_INVALID_STEP
                     end if
                  else
                     status = MC_INVALID_STEP
                  end if
               end if
            end if
         end if
      end if
      if (status /= MC_OK) T = 0.0_wp
      if (present(status_out)) then
         status_out = status
      else if (status /= MC_OK) then
         error stop 'grain_T_from_U: analytic enthalpy inversion failed'
      end if
   end function grain_T_from_U


   subroutine grain_radiative_equilibrium(g, Teq, p_abs, tau_relax, status)
      ! The prescribed field includes both continuous and sampled photons.
      ! The 0.5 K floor is a numerical limit, not a physical equilibrium.
      type(mc_grain_t), intent(in) :: g
      real(wp), intent(out) :: Teq, p_abs, tau_relax
      integer, intent(out) :: status
      real(wp) :: a_cm, integ, lo, hi, mid, p_lo, p_hi, p_mid
      real(wp) :: delta, d_pemit, p_plus, p_minus, u_eq, c_eq
      integer :: i, iter

      status = MC_INVALID_STEP
      Teq = 0.0_wp
      p_abs = 0.0_wp
      tau_relax = 0.0_wp
      if (g%a_um <= 0.0_wp .or. .not. ieee_is_finite(g%a_um)) return
      if (.not. allocated(g%lam_grid) .or. .not. allocated(g%Q_grid) .or. &
          .not. allocated(g%u_lam_grid) .or. .not. allocated(g%T_grid) .or. &
          .not. allocated(g%QPL_grid)) return
      if (g%NLAM < 2 .or. g%NT < 2) return
      if (size(g%lam_grid) < g%NLAM .or. size(g%Q_grid) < g%NLAM .or. &
          size(g%u_lam_grid) < g%NLAM .or. size(g%T_grid) < g%NT .or. &
          size(g%QPL_grid) < g%NT) return
      if (.not. all(ieee_is_finite(g%lam_grid(1:g%NLAM))) .or. &
          .not. all(ieee_is_finite(g%Q_grid(1:g%NLAM))) .or. &
          .not. all(ieee_is_finite(g%u_lam_grid(1:g%NLAM)))) return
      if (any(g%lam_grid(2:g%NLAM) <= g%lam_grid(1:g%NLAM-1)) .or. &
          any(g%Q_grid(1:g%NLAM) < 0.0_wp) .or. &
          any(g%u_lam_grid(1:g%NLAM) < 0.0_wp)) return
      a_cm = g%a_um * 1.0e-4_wp
      integ = 0.0_wp
      do i = 1, g%NLAM-1
         integ = integ + 0.5_wp * (g%lam_grid(i+1)-g%lam_grid(i)) * &
            (g%Q_grid(i)*g%u_lam_grid(i) + g%Q_grid(i+1)*g%u_lam_grid(i+1))
      end do
      p_abs = pi*a_cm**2*c_cgs*integ
      if (.not. ieee_is_finite(p_abs) .or. p_abs < 0.0_wp) return
      lo = 0.5_wp
      p_lo = 4.0_wp*pi*a_cm**2*sigma_sb_cgs*grain_get_QPL(g,lo)*lo**4
      if (.not. ieee_is_finite(p_lo) .or. p_lo <= 0.0_wp) return
      if (p_abs <= p_lo) then
         Teq = lo
      else
         hi = max(g%T_grid(g%NT),2.0_wp*lo)
         do iter = 1, 64
            p_hi = 4.0_wp*pi*a_cm**2*sigma_sb_cgs*grain_get_QPL(g,hi)*hi**4
            if (.not. ieee_is_finite(p_hi)) return
            if (p_hi >= p_abs) exit
            if (hi > huge(hi)/2.0_wp) return
            hi = 2.0_wp*hi
         end do
         if (p_hi < p_abs) return
         do iter = 1, 100
            mid = 0.5_wp*(lo+hi)
            p_mid = 4.0_wp*pi*a_cm**2*sigma_sb_cgs*grain_get_QPL(g,mid)*mid**4
            if (.not. ieee_is_finite(p_mid)) return
            if (p_mid < p_abs) then
               lo = mid
            else
               hi = mid
            end if
            if (hi-lo <= 1.0e-12_wp*max(hi,1.0_wp)) exit
         end do
         Teq = 0.5_wp*(lo+hi)
      end if
      u_eq = U_of_T(Teq,g%a_um,g%comp)
      c_eq = grain_get_C(g,Teq)
      delta = max(1.0e-4_wp*Teq,1.0e-5_wp)
      p_plus = 4.0_wp*pi*a_cm**2*sigma_sb_cgs*grain_get_QPL(g,Teq+delta)*(Teq+delta)**4
      p_minus = 4.0_wp*pi*a_cm**2*sigma_sb_cgs*grain_get_QPL(g,max(0.5_wp,Teq-delta)) * &
                max(0.5_wp,Teq-delta)**4
      d_pemit = (p_plus-p_minus)/(Teq+delta-max(0.5_wp,Teq-delta))
      if (.not. ieee_is_finite(u_eq) .or. .not. ieee_is_finite(c_eq) .or. &
          .not. ieee_is_finite(d_pemit)) return
      if (u_eq <= 0.0_wp .or. c_eq <= 0.0_wp .or. d_pemit <= 0.0_wp) return
      ! tau_relax = max(U(Teq)/Pabs, C(Teq)/(dPemit/dT)).
      ! For Pabs=0 only the local cooling scale is defined. This scale
      ! guides a finite warm-up; it cannot prove stationarity.
      tau_relax = c_eq/d_pemit
      if (p_abs > 0.0_wp) tau_relax = max(tau_relax,u_eq/p_abs)
      if (.not. ieee_is_finite(tau_relax) .or. tau_relax <= 0.0_wp) return
      status = MC_OK
   end subroutine grain_radiative_equilibrium


   subroutine mc_equilibrate(g,rng,T_initial,T_record,warm_time,warm_events,status, &
                              duration_factor,expected_events,max_events,chunk_steps,max_steps_total)
      ! Evolve before any recording begins. Only the completed state and RNG
      ! stream are returned; failed runs leave the caller's RNG unchanged.
      ! The recording boundary is fixed before sampling photons. For a constant
      ! temperature-independent Poisson rate, discarding the sampled wait that
      ! crosses this boundary and drawing a fresh wait in the recording engine
      ! gives the correct independent residual by memorylessness.
      type(mc_grain_t), intent(in) :: g
      type(rng_t), intent(inout) :: rng
      real(wp), intent(in) :: T_initial
      real(wp), intent(out) :: T_record,warm_time
      integer, intent(out) :: warm_events,status
      real(wp), intent(in), optional :: duration_factor
      integer, intent(in), optional :: expected_events,max_events,chunk_steps,max_steps_total
      type(rng_t) :: trial_rng
      real(wp) :: Teq,p_abs,tau,target_time,factor,elapsed,dt_wait,dt_left,dt_act
      real(wp) :: tcur,tnew,dts,tloc,slope,lambda,e_before,e_after
      integer :: event_expected,event_max,step_max,chunk_max,nsteps,nev,chunk_count
      logical :: event_due,at_boundary

      T_record = 0.0_wp
      warm_time = 0.0_wp
      warm_events = 0
      status = MC_INVALID_STEP
      if (.not. ieee_is_finite(T_initial) .or. T_initial < 0.5_wp) return
      call grain_radiative_equilibrium(g,Teq,p_abs,tau,status)
      if (status /= MC_OK) return
      ! Ten relaxation times is a default duration, not a convergence test.
      factor = 10.0_wp
      if (present(duration_factor)) factor = duration_factor
      event_expected = 100
      if (present(expected_events)) event_expected = expected_events
      event_max = 1000000
      if (present(max_events)) event_max = max_events
      step_max = DEFAULT_STEP_BUDGET
      if (present(max_steps_total)) step_max = max_steps_total
      chunk_max = 300
      if (present(chunk_steps)) chunk_max = chunk_steps
      status = MC_STEP_BUDGET
      if (event_expected < 0 .or. event_max < event_expected .or. step_max <= 0 .or. chunk_max <= 0) return
      status = MC_INVALID_STEP
      if (.not. ieee_is_finite(factor) .or. factor <= 0.0_wp) return
      target_time = factor*tau
      if (.not. ieee_is_finite(target_time) .or. target_time <= 0.0_wp) return
      if (.not. ieee_is_finite(g%rate_event) .or. g%rate_event < 0.0_wp) then
         status = MC_INVALID_WAIT
         return
      end if
      ! event_expected sets the mean count at the fixed recording boundary.
      if (g%rate_event > 0.0_wp) then
         target_time = max(target_time,real(event_expected,wp)/g%rate_event)
         if (.not. ieee_is_finite(target_time) .or. target_time <= 0.0_wp) then
            status = MC_INVALID_WAIT
            return
         end if
      end if
      trial_rng = rng
      tcur = T_initial
      elapsed = 0.0_wp
      nev = 0
      nsteps = 0
      status = MC_OK
      do while (elapsed < target_time)
         if (g%rate_event > 0.0_wp) then
            dt_wait = rng_exp(trial_rng,g%rate_event)
         else
            dt_wait = target_time-elapsed
         end if
         if (ieee_is_nan(dt_wait) .or. dt_wait <= 0.0_wp) then
            status = MC_INVALID_WAIT
            exit
         end if
         dt_act = min(dt_wait,target_time-elapsed)
         at_boundary = dt_wait >= target_time-elapsed
         event_due = .not. at_boundary .and. g%rate_event > 0.0_wp
         if (event_due .and. nev >= event_max) then
            status = MC_STEP_BUDGET
            exit
         end if
         dt_left = dt_act
         do while (dt_left > 0.0_wp)
            ! Continue the same waiting interval after each bounded chunk.
            do chunk_count = 1, chunk_max
               if (nsteps >= step_max) then
                  status = MC_STEP_BUDGET
                  exit
               end if
               call step_advance(g,tcur,dt_left,g%a_um*1.0e-4_wp,tnew,dts,tloc,slope)
               if (.not. ieee_is_finite(dts) .or. .not. ieee_is_finite(tnew) .or. &
                   .not. ieee_is_finite(tloc) .or. .not. ieee_is_finite(slope) .or. &
                   dts <= 0.0_wp .or. dts > dt_left .or. dt_left-dts >= dt_left .or. &
                   tnew < 0.5_wp) then
                  status = MC_INVALID_STEP
                  exit
               end if
               nsteps = nsteps+1
               tcur = tnew
               dt_left = dt_left-dts
               if (dt_left <= 0.0_wp) exit
            end do
            if (status /= MC_OK) exit
         end do
         if (status /= MC_OK) exit
         if (.not. ieee_is_finite(elapsed+dt_act) .or. elapsed+dt_act <= elapsed) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (at_boundary) then
            elapsed = target_time
            exit
         end if
         elapsed = elapsed+dt_act
         lambda = sample_photon_lam(g,trial_rng)
         if (.not. ieee_is_finite(lambda) .or. lambda <= 0.0_wp) then
            status = MC_INVALID_STEP
            exit
         end if
         e_before = U_of_T(tcur,g%a_um,g%comp)
         e_after = e_before + hc_erg_um/lambda
         tnew = grain_T_from_U(g,e_after,status)
         if (status /= MC_OK) exit
         if (tnew < tcur) then
            if (e_after-e_before > 64.0_wp*epsilon(e_before)*abs(e_before)) then
               status = MC_INVALID_STEP
               exit
            end if
            tnew = tcur
         end if
         tcur = tnew
         nev = nev+1
      end do
      if (status /= MC_OK) return
      T_record = tcur
      warm_time = elapsed
      warm_events = nev
      rng = trial_rng
   end subroutine mc_equilibrate


   ! =========================================================================
   ! MC engine: integrate cooling, sample events, accumulate P(T)
   ! =========================================================================

   subroutine mc_run_engine(g, rng, T_init, N_events, t_max, &
                            NHIST, T_edges, dP_dT, dP_dlnT, &
                            t_total, e_abs_tot, e_emit_tot, &
                            time_rec_out, temp_pre_rec_out, temp_rec_out, &
                            dtime_rec_out, n_rec_out, status_out, chunk_steps, max_steps_total, omitted_time_out, &
                            floor_time_out, out_of_grid_time_out, sample_count_out, stored_energy_out)
      ! Inputs:
      !   g          : grain state (filled by grain_setup_*)
      !   rng        : RNG state (thread-local)
      !   T_init     : initial T [K]
      !   N_events   : number of stochastic absorption events to simulate
      !   t_max      : maximum simulated time [s]; 0 -> unlimited
      !   NHIST      : number of histogram bins for P(T)
      !   T_edges(NHIST+1)   : log-spaced edges (filled here from 0.5 K to 5000 K)
      ! Outputs:
      !   dP_dT, dP_dlnT(NHIST) : probability density and dP/dlnT
      !   t_total              : total simulated time [s]
      !   e_abs_tot, e_emit_tot: integrated absorbed/emitted energy [erg]
      ! Optional outputs (for trajectory plot, e.g. PIIM Fig 24.5):
      !   time_rec_out(:)     : time at each event boundary [s]
      !   temp_pre_rec_out(:) : T just BEFORE the photon jump (= cool_end) [K]
      !   temp_rec_out(:)     : T just AFTER the photon jump [K]
      !   dtime_rec_out(:)    : cooling segment dt that preceded the record [s]
      !   n_rec_out           : number of records actually written (<= N_events+1)
      type(mc_grain_t), intent(in)    :: g
      type(rng_t),      intent(inout) :: rng
      real(wp),         intent(in)    :: T_init, t_max
      integer,          intent(in)    :: N_events, NHIST
      real(wp),         intent(out)   :: T_edges(NHIST+1)
      real(wp),         intent(out)   :: dP_dT(NHIST), dP_dlnT(NHIST)
      real(wp),         intent(out)   :: t_total, e_abs_tot, e_emit_tot
      real(wp),         intent(out), optional :: time_rec_out(:), temp_pre_rec_out(:), &
                                                 temp_rec_out(:), dtime_rec_out(:)
      integer,          intent(out), optional :: n_rec_out, status_out
      integer,          intent(in), optional :: chunk_steps, max_steps_total
      real(wp),         intent(out), optional :: omitted_time_out, floor_time_out, out_of_grid_time_out
      integer(kind=8),  intent(out), optional :: sample_count_out
      real(wp),         intent(out), optional :: stored_energy_out

      real(wp), parameter :: T_HMIN = 0.5_wp, T_HMAX = 5000.0_wp
      real(wp) :: hist_weight(NHIST)
      real(wp) :: Tcur, time_cur, dt_event, dt_act
      real(wp) :: lam_sample, U_before, U_after, T_after, dE_phot
      real(wp) :: T_cool_end, e_emit_step, lT_min, dlnT_bin, dt_done
      real(wp) :: floor_time, out_of_grid_time
      integer  :: ev, j, n_rec, n_rec_cap, status, steps_used, step_budget, chunk_size
      integer(kind=8) :: sample_count
      logical  :: do_record

      lT_min   = log(T_HMIN)
      dlnT_bin = (log(T_HMAX) - lT_min) / real(NHIST, wp)
      do j = 1, NHIST+1
         T_edges(j) = exp(lT_min + (j-1)*dlnT_bin)
      end do
      hist_weight = 0.0_wp
      dt_act = 0.0_wp
      dt_done = 0.0_wp
      floor_time = 0.0_wp
      out_of_grid_time = 0.0_wp
      sample_count = 0_8

      status = MC_OK
      if (.not. ieee_is_finite(T_init) .or. T_init < 0.5_wp) status = MC_INVALID_STEP
      steps_used = 0
      step_budget = DEFAULT_STEP_BUDGET
      if (present(max_steps_total)) step_budget = max_steps_total
      chunk_size = 4000
      if (present(chunk_steps)) chunk_size = chunk_steps
      Tcur     = T_init
      time_cur = 0.0_wp
      e_abs_tot  = 0.0_wp
      e_emit_tot = 0.0_wp
      if (present(omitted_time_out)) omitted_time_out = 0.0_wp

      do_record = present(time_rec_out) .and. present(temp_pre_rec_out) .and. &
                  present(temp_rec_out) .and. present(dtime_rec_out) .and. &
                  present(n_rec_out)
      n_rec = 0
      n_rec_cap = 0
      if (do_record) then
         n_rec_cap = min(size(time_rec_out), size(temp_pre_rec_out), &
                         size(temp_rec_out), size(dtime_rec_out))
         if (n_rec_cap > 0) then
            n_rec = 1
            time_rec_out(1)      = 0.0_wp
            temp_pre_rec_out(1)  = T_init
            temp_rec_out(1)      = T_init
            dtime_rec_out(1)     = 0.0_wp
         end if
      end if

      do ev = 1, N_events
         if (status /= MC_OK) exit
         dt_act = 0.0_wp
         dt_done = 0.0_wp
         if (.not. ieee_is_finite(g%rate_event) .or. g%rate_event < 0.0_wp .or. &
             (g%rate_event <= 0.0_wp .and. t_max <= 0.0_wp)) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (g%rate_event <= 0.0_wp) then
            dt_event = huge(1.0_wp)
         else
            dt_event = rng_exp(rng, g%rate_event)
         end if
         if (ieee_is_nan(dt_event) .or. dt_event <= 0.0_wp .or. &
             (.not. ieee_is_finite(dt_event) .and. t_max <= 0.0_wp)) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (t_max > 0.0_wp .and. dt_event >= t_max - time_cur) then
            dt_act = t_max - time_cur
            steps_used = 0
            traj_seg_start_time = time_cur
            call cool_segment(g, Tcur, dt_act, T_cool_end, e_emit_step, &
                              hist_weight, lT_min, dlnT_bin, dt_done, status, &
                              steps_used, step_budget, chunk_size, floor_time, &
                              out_of_grid_time, sample_count)
            Tcur = T_cool_end
            if (dt_done > 0.0_wp .and. time_cur + dt_done <= time_cur) status = MC_INVALID_STEP
            time_cur = time_cur + dt_done
            e_emit_tot = e_emit_tot + e_emit_step
            if (do_record .and. n_rec < n_rec_cap) then
               n_rec = n_rec + 1
               time_rec_out(n_rec)      = time_cur
               temp_pre_rec_out(n_rec)  = Tcur    ! cooling-only, no jump beyond t_max
               temp_rec_out(n_rec)      = Tcur
               dtime_rec_out(n_rec)     = dt_done
            end if
            exit
         end if
         dt_act = dt_event
         steps_used = 0
         traj_seg_start_time = time_cur
         call cool_segment(g, Tcur, dt_act, T_cool_end, e_emit_step, &
                           hist_weight, lT_min, dlnT_bin, dt_done, status, &
                              steps_used, step_budget, chunk_size, floor_time, &
                              out_of_grid_time, sample_count)
         Tcur = T_cool_end
         if (dt_done > 0.0_wp .and. time_cur + dt_done <= time_cur) status = MC_INVALID_STEP
         time_cur = time_cur + dt_done
         e_emit_tot = e_emit_tot + e_emit_step

         if (status /= MC_OK) exit
         lam_sample = sample_photon_lam(g, rng)
         dE_phot    = hc_erg_um / lam_sample
         U_before   = U_of_T(Tcur, g%a_um, g%comp)
         U_after    = U_before + dE_phot
         T_after    = grain_T_from_U(g, U_after, status)
         if (status /= MC_OK) exit
         if (T_after < Tcur) then
            if (U_after - U_before > 64.0_wp * epsilon(U_before) * abs(U_before)) then
               status = MC_INVALID_STEP
               exit
            end if
            T_after = Tcur
         end if
         e_abs_tot  = e_abs_tot + dE_phot
         Tcur       = T_after
         ! Force-push the post-jump T to capture the vertical event jump.
         call traj_push(time_cur, Tcur, .true.)

         if (do_record .and. n_rec < n_rec_cap) then
            n_rec = n_rec + 1
            time_rec_out(n_rec)      = time_cur
            temp_pre_rec_out(n_rec)  = T_cool_end   ! T just before photon jump
            temp_rec_out(n_rec)      = Tcur         ! T just after photon jump
            dtime_rec_out(n_rec)     = dt_done
         end if
      end do

      if (status == MC_OK .and. traj_enabled .and. traj_overflow) status = MC_SAMPLE_CAPACITY
      if (present(status_out)) status_out = status
      if (present(floor_time_out)) floor_time_out = floor_time
      if (present(out_of_grid_time_out)) out_of_grid_time_out = out_of_grid_time
      if (present(sample_count_out)) sample_count_out = sample_count
      if (present(omitted_time_out)) omitted_time_out = merge(max(0.0_wp, dt_act - dt_done), &
                                                             0.0_wp, status /= MC_OK .and. ev <= N_events)
      if (do_record .and. present(n_rec_out)) n_rec_out = n_rec
      t_total = time_cur
      if (present(stored_energy_out)) then
         stored_energy_out = 0.0_wp
         if (ieee_is_finite(T_init) .and. ieee_is_finite(Tcur) .and. &
             T_init >= 0.5_wp .and. Tcur >= 0.5_wp) &
            stored_energy_out = U_of_T(Tcur, g%a_um, g%comp) - U_of_T(T_init, g%a_um, g%comp)
      end if

      dP_dT   = 0.0_wp
      dP_dlnT = 0.0_wp
      if (t_total > 0.0_wp .and. status == MC_OK) then
         do j = 1, NHIST
            dP_dT(j)   = hist_weight(j) / (t_total * (T_edges(j+1) - T_edges(j)))
            dP_dlnT(j) = hist_weight(j) / (t_total * dlnT_bin)
         end do
      end if
      if (status /= MC_OK .and. .not. present(status_out)) &
         error stop 'mc_run_engine: incomplete cooling interval'
   end subroutine mc_run_engine


   subroutine cool_segment(g, T_in, dt_total, T_out, e_emit, &
                           hist_weight, lT_min, dlnT_bin, dt_done, status, &
                           steps_used, step_budget, chunk_size, floor_time, &
                           out_of_grid_time, sample_count)
      type(mc_grain_t), intent(in) :: g
      real(wp), intent(in) :: T_in, dt_total, lT_min, dlnT_bin
      real(wp), intent(out) :: T_out, e_emit, dt_done
      real(wp), intent(inout) :: hist_weight(:)
      integer, intent(inout) :: steps_used, status
      integer, intent(in) :: step_budget, chunk_size
      real(wp), intent(inout) :: floor_time, out_of_grid_time
      integer(kind=8), intent(inout) :: sample_count
      real(wp) :: floor_step, out_step
      real(wp) :: T0, T1, dts, T_eq_loc, slope, a_cm, e_step, next_time
      logical :: valid_step
      integer :: n_chunk

      a_cm = g%a_um * 1.0e-4_wp
      T0 = T_in
      dt_done = 0.0_wp
      e_emit = 0.0_wp
      T_out = T_in
      if (dt_total <= 0.0_wp) return
      call traj_push(traj_seg_start_time, T0, .true.)
      if (chunk_size <= 0 .or. step_budget <= 0) status = MC_STEP_BUDGET
      do while (dt_done < dt_total .and. status == MC_OK)
         n_chunk = 0
         do while (dt_done < dt_total .and. n_chunk < chunk_size)
            if (steps_used >= step_budget) then
               status = MC_STEP_BUDGET
               exit
            end if
            call step_advance(g, T0, dt_total - dt_done, a_cm, &
                              T1, dts, T_eq_loc, slope)
            next_time = dt_done + dts
            if (.not. ieee_is_finite(dts) .or. .not. ieee_is_finite(T1) .or. &
                .not. ieee_is_finite(T_eq_loc) .or. .not. ieee_is_finite(slope) .or. &
                .not. ieee_is_finite(next_time) .or. dts <= 0.0_wp .or. &
                next_time <= dt_done) then
               status = MC_INVALID_STEP
               exit
            end if
            call accumulate_step(g, T0, T1, T_eq_loc, slope, dts, &
                                 hist_weight, lT_min, dlnT_bin, a_cm, e_emit, &
                                 e_step, floor_step, out_step, valid_step)
            if (.not. valid_step) then
               status = MC_INVALID_STEP
               exit
            end if
            e_emit = e_emit + e_step
            floor_time = floor_time + floor_step
            out_of_grid_time = out_of_grid_time + out_step
            sample_count = sample_count + 8_8
            T0 = T1
            dt_done = next_time
            steps_used = steps_used + 1
            n_chunk = n_chunk + 1
            call traj_push(traj_seg_start_time + dt_done, T0, .false.)
            if (out_step > 0.0_wp) then
               status = MC_HISTOGRAM_RANGE
               exit
            end if
         end do
      end do
      call traj_push(traj_seg_start_time + dt_done, T0, .true.)
      T_out = T0
   end subroutine cool_segment


   subroutine accumulate_step(g, T0, T1, T_eq_loc, slope, dts, &
                              hist_weight, lT_min, dlnT_bin, a_cm, e_emit_prior, &
                              e_emit_out, floor_step, out_step, valid_out)
      ! Sub-sample the cooling step in time, assigning histogram weight to
      ! each bin and accumulating the emission integral over the step.
      ! Returns the step's emission so the caller can apply the
      ! energy-conservation correction.
      type(mc_grain_t), intent(in)    :: g
      real(wp),         intent(in)    :: T0, T1, T_eq_loc, slope, dts
      real(wp),         intent(inout) :: hist_weight(:)
      real(wp),         intent(in)    :: lT_min, dlnT_bin, a_cm
      real(wp),         intent(in)    :: e_emit_prior
      real(wp),         intent(out)   :: e_emit_out, floor_step, out_step
      logical,          intent(out)   :: valid_out
      integer,  parameter :: N_SUB = 8
      real(wp) :: dts_sub, time_sub, Tval_sub
      integer  :: i, ibin(N_SUB)
      dts_sub = dts / real(N_SUB, wp)
      e_emit_out = 0.0_wp
      floor_step = 0.0_wp
      out_step = 0.0_wp
      valid_out = .false.
      do i = 1, N_SUB
         time_sub = (real(i,wp) - 0.5_wp) * dts_sub
         if (abs(slope) > 0.0_wp .and. abs(T0 - T_eq_loc) > 0.0_wp) then
            Tval_sub = T_eq_loc + (T0 - T_eq_loc) * exp(slope * time_sub)
         else
            Tval_sub = T0 + (T1 - T0) * (time_sub / dts)
         end if
         if (.not. ieee_is_finite(Tval_sub)) return
         if (Tval_sub < 0.5_wp) Tval_sub = 0.5_wp
         if (Tval_sub <= 0.5_wp) floor_step = floor_step + dts_sub
         ibin(i) = 0
         if (Tval_sub < 0.5_wp .or. Tval_sub > 5000.0_wp) then
            out_step = out_step + dts_sub
         else
            ibin(i) = max(1, min(size(hist_weight), &
                       1 + floor((log(Tval_sub) - lT_min) / dlnT_bin)))
            if (ibin(i) < 1) out_step = out_step + dts_sub
         end if
         e_emit_out = e_emit_out + 4.0_wp * pi * a_cm**2 * grain_get_QPL(g, Tval_sub) &
                      * sigma_sb_cgs * Tval_sub**4 * dts_sub
         if (.not. ieee_is_finite(e_emit_out)) return
      end do
      if (.not. ieee_is_finite(e_emit_prior + e_emit_out)) return
      do i = 1, N_SUB
         if (ibin(i) >= 1) hist_weight(ibin(i)) = hist_weight(ibin(i)) + dts_sub
      end do
      valid_out = .true.
   end subroutine accumulate_step


   function ode_rhs(g, Tval, a_cm) result(dTdt)
      ! dT/dt = (4 pi a^2 / C_total) * (H - <Q>_T sigma T^4)
      type(mc_grain_t), intent(in) :: g
      real(wp),         intent(in) :: Tval, a_cm
      real(wp) :: dTdt, Cval
      Cval = grain_get_C(g, Tval)
      if (Cval <= 0.0_wp) then
         dTdt = 0.0_wp
         return
      end if
      dTdt = (4.0_wp * pi * a_cm**2 / Cval) * &
             (g%H_cont - grain_get_QPL(g, Tval) * sigma_sb_cgs * Tval**4)
   end function ode_rhs


   function sample_photon_lam(g, rng) result(lam)
      ! Invert the integral of the nonnegative linear photon density.
      ! Binary search selects a positive-mass interval even at a plateau;
      ! no interpolation ever crosses a gap in the illuminating field.
      type(mc_grain_t), intent(in)    :: g
      type(rng_t),      intent(inout) :: rng
      real(wp) :: lam, u, mass, width, p0, p1, disc, distance, scale
      integer  :: lo, hi, mid
      u = rng_uniform(rng)
      u = max(nearest(0.0_wp,1.0_wp), min(nearest(1.0_wp,-1.0_wp),u))
      lo = 1; hi = g%NCDF
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (g%cdf_F(mid) <= u) then; lo = mid; else; hi = mid; end if
      end do
      width = g%cdf_lam(hi)-g%cdf_lam(lo)
      mass = (u-g%cdf_F(lo))/(g%cdf_F(hi)-g%cdf_F(lo))
      p0 = g%cdf_density(lo)
      p1 = g%cdf_density(hi)
      scale = max(p0,p1)
      p0 = p0/scale
      p1 = p1/scale
      ! Scale by interval mass before the quadratic inversion, avoiding
      ! underflow for extremely weak but nonzero spectral intervals.
      if (mass <= 0.0_wp) then
         ! An exactly represented plateau probability selects the first
         ! point of the next positive interval. Its density can be zero;
         ! evaluate the endpoint directly instead of forming zero/zero.
         distance = 0.0_wp
      else
         disc = sqrt(max(0.0_wp, p0*p0 + (p1-p0)*(p0+p1)*mass))
         distance = (p0+p1)*mass/(p0+disc)
      end if
      distance = max(nearest(0.0_wp,1.0_wp), min(nearest(1.0_wp,-1.0_wp),distance))
      lam = g%cdf_lam(lo)+width*distance
      if (lam <= g%cdf_lam(lo) .and. p0 == 0.0_wp) lam = nearest(g%cdf_lam(lo),1.0_wp)
      if (lam >= g%cdf_lam(hi) .and. p1 == 0.0_wp) lam = nearest(g%cdf_lam(hi),-1.0_wp)
   end function sample_photon_lam


   ! ==========================================================================
   ! Shared helpers for adaptive-grid engines.
   ! step_advance + sub_sample_step factor out the trajectory math AND the
   ! step-size control of cool_segment so that the binning policy (fixed-grid
   ! / probe / replay / buffer) is the only thing that differs between the
   ! fixed-grid and adaptive-grid engines; the integrated trajectory is
   ! identical.
   ! ==========================================================================

   subroutine step_advance(g, T0, dt_left, a_cm, T1, dts, T_eq_loc, slope)
      ! One internal cooling step for the adaptive-grid engines.  Same
      ! trajectory math as the inline block in cool_segment() (local
      ! linearization, exponential relaxation toward T_eq) AND the same step
      ! tolerance: both the cooling and heating branches cap the fractional
      ! temperature change per step at DT_FRAC_MAX = 0.05, identical to
      ! cool_segment.  The adaptive-grid engines therefore integrate the
      ! same trajectory as the fixed-grid engine; only the histogram binning
      ! policy differs.
      type(mc_grain_t), intent(in)  :: g
      real(wp),         intent(in)  :: T0, dt_left, a_cm
      real(wp),         intent(out) :: T1, dts, T_eq_loc, slope
      real(wp) :: dTdt0, dTdt1, T_probe, half_life, rel_excursion
      real(wp), parameter :: T_FLOOR    = 0.5_wp
      real(wp), parameter :: STEP_FRAC  = 1.0_wp     ! half-life cap (coarse)
      real(wp), parameter :: DT_FRAC_MAX = 0.05_wp   ! max fractional dT per step
                                                     ! (bounds linearization error
                                                     !  of the super-linear cooling)
      real(wp), parameter :: TINY_RATE  = 1.0e-300_wp
      real(wp), parameter :: CONVERGED  = 1.0e-3_wp
      real(wp), parameter :: LN2        = 0.6931471805599453_wp

      dTdt0   = ode_rhs(g, T0, a_cm)
      ! The imposed 0.5 K numerical floor approximates the unresolved cold
      ! tail only while the physical derivative points below the floor.
      ! Include the entire remaining residence and emission at this floor.
      ! Use this cold-tail approximation only when its emission is negligible
      ! in the wavelength range being studied; it does not predict cold P(T).
      ! The caller can inspect the separate floor_time_out diagnostic.
      if (T0 <= T_FLOOR .and. dTdt0 <= 0.0_wp) then
         dts = dt_left
         T1 = T_FLOOR
         T_eq_loc = T_FLOOR
         slope = 0.0_wp
         return
      end if
      T_probe = T0 * 1.01_wp
      if (T_probe <= T0) T_probe = T0 + 1.0e-4_wp
      dTdt1   = ode_rhs(g, T_probe, a_cm)
      slope   = (dTdt1 - dTdt0) / (T_probe - T0)

      if (slope < 0.0_wp) then
         T_eq_loc      = T0 - dTdt0 / slope
         rel_excursion = abs(T0 - T_eq_loc) / max(T0, T_FLOOR)
         half_life     = LN2 / abs(slope)
         if (rel_excursion < CONVERGED) then
            dts = dt_left
            T1  = T_eq_loc
         else
            dts = min(dt_left, STEP_FRAC * half_life, &
                      DT_FRAC_MAX * T0 / max(abs(dTdt0), TINY_RATE))
            T1  = T_eq_loc + (T0 - T_eq_loc) * exp(slope * dts)
         end if
      else if (abs(dTdt0) > 0.0_wp) then
         dts      = min(dt_left, DT_FRAC_MAX * T0 / abs(dTdt0))
         T1       = T0 + dts * dTdt0
         T_eq_loc = T0
      else
         dts      = dt_left
         T1       = T0
         T_eq_loc = T0
      end if
      if (T1 < T_FLOOR) T1 = T_FLOOR
   end subroutine step_advance


   subroutine sub_sample_step(T0, T1, T_eq_loc, slope, dts, N_SUB, T_sub, dt_sub)
      ! Generate N_SUB midpoint sub-samples of the exponential trajectory
      ! T(t) = T_eq_loc + (T0 - T_eq_loc) * exp(slope * t) on [0, dts].
      ! Falls back to linear interpolation if the slope/excursion is zero.
      real(wp), intent(in)  :: T0, T1, T_eq_loc, slope, dts
      integer,  intent(in)  :: N_SUB
      real(wp), intent(out) :: T_sub(N_SUB), dt_sub
      integer  :: i
      real(wp) :: time_sub
      dt_sub = dts / real(N_SUB, wp)
      do i = 1, N_SUB
         time_sub = (real(i,wp) - 0.5_wp) * dt_sub
         if (abs(slope) > 0.0_wp .and. abs(T0 - T_eq_loc) > 0.0_wp) then
            T_sub(i) = T_eq_loc + (T0 - T_eq_loc) * exp(slope * time_sub)
         else
            T_sub(i) = T0 + (T1 - T0) * (time_sub / dts)
         end if
         if (T_sub(i) < 0.5_wp) T_sub(i) = 0.5_wp
      end do
   end subroutine sub_sample_step


   function e_emit_substep(g, T_sub, N_SUB, dt_sub, a_cm) result(e_emit)
      ! Trapezoid (midpoint) integral of 4 pi a^2 <Q>_T sigma T^4 over sub-step.
      type(mc_grain_t), intent(in) :: g
      integer,          intent(in) :: N_SUB
      real(wp),         intent(in) :: T_sub(N_SUB), dt_sub, a_cm
      real(wp) :: e_emit
      integer  :: i
      e_emit = 0.0_wp
      do i = 1, N_SUB
         e_emit = e_emit + 4.0_wp * pi * a_cm**2 * grain_get_QPL(g, T_sub(i)) &
                  * sigma_sb_cgs * T_sub(i)**4 * dt_sub
      end do
   end function e_emit_substep


   subroutine build_adaptive_grid(T_min_obs, T_max_obs, NHIST, margin, &
                                  T_edges, is_log, lT_or_T_lo, dlnT_or_dT)
      ! Build NHIST bin edges over [T_min_obs*(1-margin), T_max_obs*(1+margin)].
      ! Spacing is log if T_hi/T_lo > LOG_LINEAR_THRESHOLD, else linear --
      ! large grains with narrow Teq distributions get linear bins;
      ! small grains with multi-decade T excursions get log bins.
      integer,  intent(in)  :: NHIST
      real(wp), intent(in)  :: T_min_obs, T_max_obs, margin
      real(wp), intent(out) :: T_edges(NHIST+1), lT_or_T_lo, dlnT_or_dT
      logical,  intent(out) :: is_log
      real(wp), parameter   :: T_FLOOR = 0.5_wp
      real(wp), parameter   :: LOG_LINEAR_THRESHOLD = 3.0_wp
      real(wp) :: T_lo, T_hi
      integer  :: j

      T_lo = max(T_min_obs * (1.0_wp - margin), T_FLOOR)
      T_hi = T_max_obs * (1.0_wp + margin)
      if (T_hi <= T_lo) T_hi = T_lo * 1.10_wp + 1.0e-3_wp
      is_log = (T_hi / T_lo > LOG_LINEAR_THRESHOLD)
      if (is_log) then
         lT_or_T_lo = log(T_lo)
         dlnT_or_dT = (log(T_hi) - lT_or_T_lo) / real(NHIST, wp)
         do j = 1, NHIST+1
            T_edges(j) = exp(lT_or_T_lo + (j-1)*dlnT_or_dT)
         end do
      else
         lT_or_T_lo = T_lo
         dlnT_or_dT = (T_hi - T_lo) / real(NHIST, wp)
         do j = 1, NHIST+1
            T_edges(j) = T_lo + (j-1)*dlnT_or_dT
         end do
      end if
   end subroutine build_adaptive_grid


   subroutine bin_sub_samples(T_sub, N_SUB, dt_sub, hist_weight, NHIST, &
                              is_log, lT_or_T_lo, dlnT_or_dT, out_step, floor_step)
      ! Bin N_SUB equal-dt sub-samples into hist_weight on the adaptive grid.
      integer,  intent(in)    :: N_SUB, NHIST
      real(wp), intent(in)    :: T_sub(N_SUB), dt_sub, lT_or_T_lo, dlnT_or_dT
      logical,  intent(in)    :: is_log
      real(wp), intent(inout) :: hist_weight(NHIST)
      real(wp), intent(out) :: out_step, floor_step
      integer  :: i, ibin
      real(wp) :: x, x_hi
      out_step = 0.0_wp
      floor_step = 0.0_wp
      x_hi = lT_or_T_lo + NHIST * dlnT_or_dT
      do i = 1, N_SUB
         if (is_log) then
            x = log(max(T_sub(i), 1.0e-30_wp))
         else
            x = T_sub(i)
         end if
         if (T_sub(i) <= 0.5_wp) floor_step = floor_step + dt_sub
         if (x < lT_or_T_lo .or. x > x_hi) then
            out_step = out_step + dt_sub
         else
            ibin = max(1, min(NHIST, 1 + floor((x - lT_or_T_lo) / dlnT_or_dT)))
            hist_weight(ibin) = hist_weight(ibin) + dt_sub
         end if
      end do
   end subroutine bin_sub_samples


   ! ==========================================================================
   ! Method A: two-pass with deterministic RNG replay.
   ! Pass 1 propagates the trajectory and tracks T_min/T_max of all sub-step
   ! samples. Pass 2 restores the RNG state and re-runs with
   ! the adaptive grid built from the observed range.  Cost ~ 2 x mc_run_engine.
   ! ==========================================================================

   subroutine mc_run_engine_2pass(g, rng, T_init, N_events, N_pass1, t_max, &
                                  NHIST, T_edges, dP_dT, dP_dlnT, &
                                  t_total, e_abs_tot, e_emit_tot, &
                                  is_log_out, T_min_obs_out, T_max_obs_out, &
                                  status_out, chunk_steps, max_steps_total, omitted_time_out, &
                                  floor_time_out, out_of_grid_time_out, sample_count_out, stored_energy_out)
      ! Inputs:
      !   N_pass1 : number of events in Pass 1 (the T_min/T_max probe).
      !             Pass 2 always uses the full N_events for the binned
      !             trajectory. If N_pass1 <= 0, defaults to N_events.
      !             If N_pass1 > N_events, silently clamped to N_events.
      ! Outputs (optional):
      !   is_log_out, T_min_obs_out, T_max_obs_out : grid diagnostics
      type(mc_grain_t), intent(in)    :: g
      type(rng_t),      intent(inout) :: rng
      real(wp),         intent(in)    :: T_init, t_max
      integer,          intent(in)    :: N_events, N_pass1, NHIST
      real(wp),         intent(out)   :: T_edges(NHIST+1)
      real(wp),         intent(out)   :: dP_dT(NHIST), dP_dlnT(NHIST)
      real(wp),         intent(out)   :: t_total, e_abs_tot, e_emit_tot
      logical,          intent(out), optional :: is_log_out
      real(wp),         intent(out), optional :: T_min_obs_out, T_max_obs_out
      integer,          intent(out), optional :: status_out
      integer,          intent(in), optional :: chunk_steps, max_steps_total
      real(wp),         intent(out), optional :: omitted_time_out, floor_time_out, out_of_grid_time_out
      integer(kind=8),  intent(out), optional :: sample_count_out
      real(wp),         intent(out), optional :: stored_energy_out

      integer,  parameter :: N_SUB     = 8
      real(wp), parameter :: MARGIN    = 0.05_wp

      integer(kind=8) :: rng_saved
      real(wp) :: hist_weight(NHIST)
      real(wp) :: Tcur, time_cur, dt_event, dt_act, a_cm, tau_seg, omitted_time
      real(wp) :: lam_sample, U_before, U_after, dE_phot
      real(wp) :: T_min_obs, T_max_obs, floor_time, out_of_grid_time, floor_step, out_step
      real(wp) :: T0_step, T1_step, dts, T_eq_loc, slope
      real(wp) :: T_sub(N_SUB), dt_sub, e_step, lT_or_T_lo, dlnT_or_dT
      logical  :: is_log
      integer  :: ev, j, nstep, N_pass1_eff, status, steps_used, step_budget, chunk_size
      integer(kind=8) :: nseg_call, ncap_call, sample_count

      a_cm = g%a_um * 1.0e-4_wp
      nseg_call = 0_8
      ncap_call = 0_8
      omitted_time = 0.0_wp
      floor_time = 0.0_wp
      out_of_grid_time = 0.0_wp
      sample_count = 0_8
      dt_act = 0.0_wp
      status = MC_OK
      if (.not. ieee_is_finite(T_init) .or. T_init < 0.5_wp) status = MC_INVALID_STEP
      steps_used = 0
      step_budget = DEFAULT_STEP_BUDGET
      if (present(max_steps_total)) step_budget = max_steps_total
      chunk_size = 300
      if (present(chunk_steps)) chunk_size = chunk_steps
      if (step_budget <= 0 .or. chunk_size <= 0) status = MC_STEP_BUDGET
      N_pass1_eff = N_pass1
      ! A shorter probe can miss later temperatures. Replay reports
      ! MC_HISTOGRAM_RANGE if any sample falls outside its grid.
      if (N_pass1_eff <= 0)         N_pass1_eff = N_events
      if (N_pass1_eff > N_events)   N_pass1_eff = N_events

      rng_saved = rng%state

      ! ---- Pass 1: probe (only first N_pass1_eff events) ----
      Tcur     = T_init
      time_cur = 0.0_wp
      T_min_obs = huge(1.0_wp)
      T_max_obs = 0.0_wp
      do ev = 1, N_pass1_eff
         if (status /= MC_OK) exit
         if (.not. ieee_is_finite(g%rate_event) .or. g%rate_event < 0.0_wp .or. &
             (g%rate_event <= 0.0_wp .and. t_max <= 0.0_wp)) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (g%rate_event <= 0.0_wp) then
            dt_event = huge(1.0_wp)
         else
            dt_event = rng_exp(rng, g%rate_event)
         end if
         if (ieee_is_nan(dt_event) .or. dt_event <= 0.0_wp .or. &
             (.not. ieee_is_finite(dt_event) .and. t_max <= 0.0_wp)) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (t_max > 0.0_wp .and. dt_event >= t_max - time_cur) then
            dt_act = t_max - time_cur
         else
            dt_act = dt_event
         end if

         steps_used = 0
         T0_step = Tcur
         tau_seg = 0.0_wp
         nstep   = 0
         do while (tau_seg < dt_act .and. status == MC_OK)
            if (steps_used >= step_budget) then
               status = MC_STEP_BUDGET
               exit
            end if
            call step_advance(g, T0_step, dt_act - tau_seg, a_cm, &
                              T1_step, dts, T_eq_loc, slope)
            if (.not. ieee_is_finite(dts) .or. .not. ieee_is_finite(T1_step) .or. &
                .not. ieee_is_finite(T_eq_loc) .or. .not. ieee_is_finite(slope) .or. &
                .not. ieee_is_finite(tau_seg + dts) .or. dts <= 0.0_wp .or. &
                tau_seg + dts <= tau_seg) then
               status = MC_INVALID_STEP
               exit
            end if
            steps_used = steps_used + 1
            nstep = nstep + 1
            if (mod(nstep, chunk_size) == 0 .and. tau_seg + dts < dt_act) &
               ncap_call = ncap_call + 1_8
            call sub_sample_step(T0_step, T1_step, T_eq_loc, slope, dts, &
                                 N_SUB, T_sub, dt_sub)
            if (.not. all(ieee_is_finite(T_sub))) then
               status = MC_INVALID_STEP
               exit
            end if
            sample_count = sample_count + int(N_SUB, kind=8)
            floor_time = floor_time + count(T_sub <= 0.5_wp) * dt_sub
            T_min_obs = min(T_min_obs, minval(T_sub))
            T_max_obs = max(T_max_obs, maxval(T_sub))
            T0_step = T1_step
            tau_seg = tau_seg + dts
         end do
         nseg_call = nseg_call + 1_8
         Tcur     = T0_step
         if (tau_seg > 0.0_wp .and. time_cur + tau_seg <= time_cur) status = MC_INVALID_STEP
         time_cur = time_cur + tau_seg
         if (status /= MC_OK) omitted_time = max(0.0_wp, dt_act - tau_seg)

         if (status /= MC_OK) exit
         if (t_max > 0.0_wp .and. time_cur >= t_max) exit

         lam_sample = sample_photon_lam(g, rng)
         dE_phot    = hc_erg_um / lam_sample
         U_before   = U_of_T(Tcur, g%a_um, g%comp)
         U_after    = U_before + dE_phot
         T1_step = grain_T_from_U(g, U_after, status)
         if (status /= MC_OK) exit
         if (T1_step < Tcur) then
            if (U_after - U_before > 64.0_wp * epsilon(U_before) * abs(U_before)) then
               status = MC_INVALID_STEP
               exit
            end if
            T1_step = Tcur
         end if
         Tcur = T1_step
      end do

      if (T_min_obs >= T_max_obs) then
         T_min_obs = max(T_init * 0.5_wp, 0.5_wp)
         T_max_obs = T_init * 2.0_wp + 1.0_wp
      end if

      call build_adaptive_grid(T_min_obs, T_max_obs, NHIST, MARGIN, &
                               T_edges, is_log, lT_or_T_lo, dlnT_or_dT)
      if (present(is_log_out))    is_log_out    = is_log
      if (present(T_min_obs_out)) T_min_obs_out = T_min_obs
      if (present(T_max_obs_out)) T_max_obs_out = T_max_obs

      if (status /= MC_OK) then
         ! The probe has no emission integral. Its failure diagnostics report
         ! elapsed probe time and omitted time; production energies are unset.
         t_total = time_cur
         if (present(stored_energy_out)) then
            stored_energy_out = 0.0_wp
            if (ieee_is_finite(T_init) .and. ieee_is_finite(Tcur) .and. &
                T_init >= 0.5_wp .and. Tcur >= 0.5_wp) &
               stored_energy_out = U_of_T(Tcur, g%a_um, g%comp) - U_of_T(T_init, g%a_um, g%comp)
         end if
         e_abs_tot = 0.0_wp
         e_emit_tot = 0.0_wp
         dP_dT = 0.0_wp
         dP_dlnT = 0.0_wp
         if (present(status_out)) status_out = status
         if (present(omitted_time_out)) omitted_time_out = omitted_time
         if (present(floor_time_out)) floor_time_out = floor_time
         if (present(out_of_grid_time_out)) out_of_grid_time_out = out_of_grid_time
         if (present(sample_count_out)) sample_count_out = sample_count
         call stepcap_tally('2pass', g%a_um, nseg_call, ncap_call)
         if (.not. present(status_out)) &
            error stop 'mc_run_engine_2pass: incomplete probe interval'
         return
      end if
      steps_used = 0
      omitted_time = 0.0_wp
      floor_time = 0.0_wp
      out_of_grid_time = 0.0_wp
      sample_count = 0_8

      ! ---- Pass 2: replay with binning ----
      rng%state   = rng_saved
      Tcur        = T_init
      time_cur    = 0.0_wp
      e_abs_tot   = 0.0_wp
      e_emit_tot  = 0.0_wp
      hist_weight = 0.0_wp

      do ev = 1, N_events
         if (status /= MC_OK) exit
         if (.not. ieee_is_finite(g%rate_event) .or. g%rate_event < 0.0_wp .or. &
             (g%rate_event <= 0.0_wp .and. t_max <= 0.0_wp)) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (g%rate_event <= 0.0_wp) then
            dt_event = huge(1.0_wp)
         else
            dt_event = rng_exp(rng, g%rate_event)
         end if
         if (ieee_is_nan(dt_event) .or. dt_event <= 0.0_wp .or. &
             (.not. ieee_is_finite(dt_event) .and. t_max <= 0.0_wp)) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (t_max > 0.0_wp .and. dt_event >= t_max - time_cur) then
            dt_act = t_max - time_cur
         else
            dt_act = dt_event
         end if

         steps_used = 0
         T0_step = Tcur
         tau_seg = 0.0_wp
         nstep   = 0
         do while (tau_seg < dt_act .and. status == MC_OK)
            if (steps_used >= step_budget) then
               status = MC_STEP_BUDGET
               exit
            end if
            call step_advance(g, T0_step, dt_act - tau_seg, a_cm, &
                              T1_step, dts, T_eq_loc, slope)
            if (.not. ieee_is_finite(dts) .or. .not. ieee_is_finite(T1_step) .or. &
                .not. ieee_is_finite(T_eq_loc) .or. .not. ieee_is_finite(slope) .or. &
                .not. ieee_is_finite(tau_seg + dts) .or. dts <= 0.0_wp .or. &
                tau_seg + dts <= tau_seg) then
               status = MC_INVALID_STEP
               exit
            end if
            steps_used = steps_used + 1
            nstep = nstep + 1
            if (mod(nstep, chunk_size) == 0 .and. tau_seg + dts < dt_act) &
               ncap_call = ncap_call + 1_8
            call sub_sample_step(T0_step, T1_step, T_eq_loc, slope, dts, &
                                 N_SUB, T_sub, dt_sub)
            if (.not. all(ieee_is_finite(T_sub))) then
               status = MC_INVALID_STEP
               exit
            end if
            e_step = e_emit_substep(g, T_sub, N_SUB, dt_sub, a_cm)
            if (.not. ieee_is_finite(e_step) .or. .not. ieee_is_finite(e_emit_tot + e_step)) then
               status = MC_INVALID_STEP
               exit
            end if
            call bin_sub_samples(T_sub, N_SUB, dt_sub, hist_weight, NHIST, &
                                 is_log, lT_or_T_lo, dlnT_or_dT, out_step, floor_step)
            out_of_grid_time = out_of_grid_time + out_step
            floor_time = floor_time + floor_step
            sample_count = sample_count + int(N_SUB, kind=8)
            e_emit_tot = e_emit_tot + e_step
            T0_step = T1_step
            tau_seg = tau_seg + dts
            if (out_step > 0.0_wp) then
               status = MC_HISTOGRAM_RANGE
               exit
            end if
         end do
         nseg_call = nseg_call + 1_8
         Tcur     = T0_step
         if (tau_seg > 0.0_wp .and. time_cur + tau_seg <= time_cur) status = MC_INVALID_STEP
         time_cur = time_cur + tau_seg
         if (status /= MC_OK) omitted_time = max(0.0_wp, dt_act - tau_seg)

         if (status /= MC_OK) exit
         if (t_max > 0.0_wp .and. time_cur >= t_max) exit

         lam_sample = sample_photon_lam(g, rng)
         dE_phot    = hc_erg_um / lam_sample
         U_before   = U_of_T(Tcur, g%a_um, g%comp)
         U_after    = U_before + dE_phot
         T1_step = grain_T_from_U(g, U_after, status)
         if (status /= MC_OK) exit
         if (T1_step < Tcur) then
            if (U_after - U_before > 64.0_wp * epsilon(U_before) * abs(U_before)) then
               status = MC_INVALID_STEP
               exit
            end if
            T1_step = Tcur
         end if
         Tcur = T1_step
         e_abs_tot = e_abs_tot + dE_phot
      end do

      t_total = time_cur
      if (present(stored_energy_out)) then
         stored_energy_out = 0.0_wp
         if (ieee_is_finite(T_init) .and. ieee_is_finite(Tcur) .and. &
             T_init >= 0.5_wp .and. Tcur >= 0.5_wp) &
            stored_energy_out = U_of_T(Tcur, g%a_um, g%comp) - U_of_T(T_init, g%a_um, g%comp)
      end if

      call stepcap_tally('2pass', g%a_um, nseg_call, ncap_call)
      if (present(status_out)) status_out = status
      if (present(omitted_time_out)) omitted_time_out = omitted_time
      if (present(floor_time_out)) floor_time_out = floor_time
      if (present(out_of_grid_time_out)) out_of_grid_time_out = out_of_grid_time
      if (present(sample_count_out)) sample_count_out = sample_count

      dP_dT   = 0.0_wp
      dP_dlnT = 0.0_wp
      if (t_total > 0.0_wp .and. status == MC_OK) then
         do j = 1, NHIST
            dP_dT(j)   = hist_weight(j) / (t_total * (T_edges(j+1) - T_edges(j)))
            dP_dlnT(j) = hist_weight(j) / (t_total * &
                         (log(T_edges(j+1)) - log(T_edges(j))))
         end do
      end if
      if (status /= MC_OK .and. .not. present(status_out)) &
         error stop 'mc_run_engine_2pass: incomplete cooling interval'
   end subroutine mc_run_engine_2pass


   ! ==========================================================================
   ! Method B: single-pass with sub-step sample buffer.
   ! Trajectory and emission are computed in one pass; sub-step (T, dt) pairs
   ! are appended to caller-bounded arrays.  T_min/T_max
   ! is tracked from all completed samples. On success the adaptive grid
   ! is built and the complete buffer is binned.
   ! ==========================================================================

   subroutine mc_run_engine_buffered(g, rng, T_init, N_events, t_max, &
                                     NHIST, max_samples, &
                                     T_edges, dP_dT, dP_dlnT, &
                                     t_total, e_abs_tot, e_emit_tot, &
                                     is_log_out, T_min_obs_out, T_max_obs_out, &
                                     n_samples_used, buffer_full_out, &
                                     status_out, chunk_steps, max_steps_total, omitted_time_out, &
                                  floor_time_out, out_of_grid_time_out, sample_count_out, stored_energy_out)
      ! Inputs:
      !   max_samples : capacity of internal sub-step (T,dt) buffer.
      !                 Insufficient capacity returns MC_SAMPLE_CAPACITY and
      !                 sample_count_out gives a lower bound. A 2-pass run
      !                 with the same RNG state gives the exact capacity.
      type(mc_grain_t), intent(in)    :: g
      type(rng_t),      intent(inout) :: rng
      real(wp),         intent(in)    :: T_init, t_max
      integer,          intent(in)    :: N_events, NHIST, max_samples
      real(wp),         intent(out)   :: T_edges(NHIST+1)
      real(wp),         intent(out)   :: dP_dT(NHIST), dP_dlnT(NHIST)
      real(wp),         intent(out)   :: t_total, e_abs_tot, e_emit_tot
      logical,          intent(out), optional :: is_log_out, buffer_full_out
      real(wp),         intent(out), optional :: T_min_obs_out, T_max_obs_out
      integer,          intent(out), optional :: status_out
      integer,          intent(in), optional :: chunk_steps, max_steps_total
      real(wp),         intent(out), optional :: omitted_time_out, floor_time_out, out_of_grid_time_out
      integer(kind=8),  intent(out), optional :: sample_count_out
      real(wp),         intent(out), optional :: stored_energy_out
      integer,          intent(out), optional :: n_samples_used

      integer,  parameter :: N_SUB     = 8
      real(wp), parameter :: MARGIN    = 0.05_wp

      real(wp) :: hist_weight(NHIST)
      real(wp), allocatable :: T_buf(:), dt_buf(:)
      real(wp) :: Tcur, time_cur, dt_event, dt_act, a_cm, tau_seg, omitted_time
      real(wp) :: lam_sample, U_before, U_after, dE_phot
      real(wp) :: T_min_obs, T_max_obs, floor_time, out_of_grid_time, floor_step, out_step
      real(wp) :: T0_step, T1_step, dts, T_eq_loc, slope
      real(wp) :: T_sub(N_SUB), dt_sub, e_step, lT_or_T_lo, dlnT_or_dT, t_buf_total
      logical  :: is_log, buf_full
      integer  :: ev, j, k, nstep, n_used, status, steps_used, step_budget, chunk_size
      integer(kind=8) :: nseg_call, ncap_call, sample_count

      a_cm = g%a_um * 1.0e-4_wp
      nseg_call = 0_8
      ncap_call = 0_8
      omitted_time = 0.0_wp
      floor_time = 0.0_wp
      out_of_grid_time = 0.0_wp
      sample_count = 0_8
      dt_act = 0.0_wp
      status = MC_OK
      if (.not. ieee_is_finite(T_init) .or. T_init < 0.5_wp) status = MC_INVALID_STEP
      steps_used = 0
      step_budget = DEFAULT_STEP_BUDGET
      if (present(max_steps_total)) step_budget = max_steps_total
      chunk_size = 300
      if (present(chunk_steps)) chunk_size = chunk_steps
      if (step_budget <= 0 .or. chunk_size <= 0) status = MC_STEP_BUDGET
      if (max_samples < 0) status = MC_SAMPLE_CAPACITY

      allocate(T_buf(max(0,max_samples)), dt_buf(max(0,max_samples)))
      n_used   = 0
      buf_full = .false.

      Tcur       = T_init
      time_cur   = 0.0_wp
      e_abs_tot  = 0.0_wp
      e_emit_tot = 0.0_wp
      T_min_obs  = huge(1.0_wp)
      T_max_obs  = 0.0_wp

      do ev = 1, N_events
         if (status /= MC_OK) exit
         if (.not. ieee_is_finite(g%rate_event) .or. g%rate_event < 0.0_wp .or. &
             (g%rate_event <= 0.0_wp .and. t_max <= 0.0_wp)) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (g%rate_event <= 0.0_wp) then
            dt_event = huge(1.0_wp)
         else
            dt_event = rng_exp(rng, g%rate_event)
         end if
         if (ieee_is_nan(dt_event) .or. dt_event <= 0.0_wp .or. &
             (.not. ieee_is_finite(dt_event) .and. t_max <= 0.0_wp)) then
            status = MC_INVALID_WAIT
            exit
         end if
         if (t_max > 0.0_wp .and. dt_event >= t_max - time_cur) then
            dt_act = t_max - time_cur
         else
            dt_act = dt_event
         end if

         steps_used = 0
         T0_step = Tcur
         tau_seg = 0.0_wp
         nstep   = 0
         do while (tau_seg < dt_act .and. status == MC_OK)
            if (steps_used >= step_budget) then
               status = MC_STEP_BUDGET
               exit
            end if
            call step_advance(g, T0_step, dt_act - tau_seg, a_cm, &
                              T1_step, dts, T_eq_loc, slope)
            if (.not. ieee_is_finite(dts) .or. .not. ieee_is_finite(T1_step) .or. &
                .not. ieee_is_finite(T_eq_loc) .or. .not. ieee_is_finite(slope) .or. &
                .not. ieee_is_finite(tau_seg + dts) .or. dts <= 0.0_wp .or. &
                tau_seg + dts <= tau_seg) then
               status = MC_INVALID_STEP
               exit
            end if
            steps_used = steps_used + 1
            nstep = nstep + 1
            if (mod(nstep, chunk_size) == 0 .and. tau_seg + dts < dt_act) &
               ncap_call = ncap_call + 1_8
            call sub_sample_step(T0_step, T1_step, T_eq_loc, slope, dts, &
                                 N_SUB, T_sub, dt_sub)
            if (.not. all(ieee_is_finite(T_sub))) then
               status = MC_INVALID_STEP
               exit
            end if
            e_step = e_emit_substep(g, T_sub, N_SUB, dt_sub, a_cm)
            if (.not. ieee_is_finite(e_step) .or. .not. ieee_is_finite(e_emit_tot + e_step)) then
               status = MC_INVALID_STEP
               exit
            end if
            if (sample_count + int(N_SUB, kind=8) > int(max_samples, kind=8)) then
               ! The next complete step would exceed capacity. Count its
               ! required samples as a lower bound without advancing time.
               sample_count = sample_count + int(N_SUB, kind=8)
               buf_full = .true.
               status = MC_SAMPLE_CAPACITY
               exit
            end if
            T_min_obs = min(T_min_obs, minval(T_sub))
            T_max_obs = max(T_max_obs, maxval(T_sub))
            floor_time = floor_time + count(T_sub <= 0.5_wp) * dt_sub
            sample_count = sample_count + int(N_SUB, kind=8)
            do k = 1, N_SUB
               n_used         = n_used + 1
               T_buf(n_used)  = T_sub(k)
               dt_buf(n_used) = dt_sub
            end do
            e_emit_tot = e_emit_tot + e_step
            T0_step = T1_step
            tau_seg = tau_seg + dts
         end do
         nseg_call = nseg_call + 1_8
         Tcur     = T0_step
         if (tau_seg > 0.0_wp .and. time_cur + tau_seg <= time_cur) status = MC_INVALID_STEP
         time_cur = time_cur + tau_seg
         if (status /= MC_OK) omitted_time = max(0.0_wp, dt_act - tau_seg)

         if (status /= MC_OK) exit
         if (t_max > 0.0_wp .and. time_cur >= t_max) exit

         lam_sample = sample_photon_lam(g, rng)
         dE_phot    = hc_erg_um / lam_sample
         U_before   = U_of_T(Tcur, g%a_um, g%comp)
         U_after    = U_before + dE_phot
         T1_step = grain_T_from_U(g, U_after, status)
         if (status /= MC_OK) exit
         if (T1_step < Tcur) then
            if (U_after - U_before > 64.0_wp * epsilon(U_before) * abs(U_before)) then
               status = MC_INVALID_STEP
               exit
            end if
            T1_step = Tcur
         end if
         Tcur = T1_step
         e_abs_tot = e_abs_tot + dE_phot
      end do

      t_total = time_cur
      if (present(stored_energy_out)) then
         stored_energy_out = 0.0_wp
         if (ieee_is_finite(T_init) .and. ieee_is_finite(Tcur) .and. &
             T_init >= 0.5_wp .and. Tcur >= 0.5_wp) &
            stored_energy_out = U_of_T(Tcur, g%a_um, g%comp) - U_of_T(T_init, g%a_um, g%comp)
      end if

      call stepcap_tally('buffered', g%a_um, nseg_call, ncap_call)
      if (present(status_out)) status_out = status
      if (present(omitted_time_out)) omitted_time_out = omitted_time
      if (present(floor_time_out)) floor_time_out = floor_time
      if (present(sample_count_out)) sample_count_out = sample_count

      if (T_min_obs >= T_max_obs) then
         T_min_obs = max(T_init * 0.5_wp, 0.5_wp)
         T_max_obs = T_init * 2.0_wp + 1.0_wp
      end if

      call build_adaptive_grid(T_min_obs, T_max_obs, NHIST, MARGIN, &
                               T_edges, is_log, lT_or_T_lo, dlnT_or_dT)
      if (present(is_log_out))     is_log_out     = is_log
      if (present(buffer_full_out)) buffer_full_out = buf_full
      if (present(T_min_obs_out))  T_min_obs_out  = T_min_obs
      if (present(T_max_obs_out))  T_max_obs_out  = T_max_obs

      hist_weight = 0.0_wp
      t_buf_total = 0.0_wp
      do k = 1, n_used
         call bin_sub_samples(T_buf(k:k), 1, dt_buf(k), hist_weight, NHIST, &
                              is_log, lT_or_T_lo, dlnT_or_dT, out_step, floor_step)
         out_of_grid_time = out_of_grid_time + out_step
         t_buf_total = t_buf_total + dt_buf(k)
      end do

      dP_dT   = 0.0_wp
      dP_dlnT = 0.0_wp
      if (present(out_of_grid_time_out)) out_of_grid_time_out = out_of_grid_time
      if (out_of_grid_time > 0.0_wp .and. status == MC_OK) status = MC_HISTOGRAM_RANGE
      if (present(status_out)) status_out = status
      if (t_buf_total > 0.0_wp .and. status == MC_OK) then
         do j = 1, NHIST
            dP_dT(j)   = hist_weight(j) / (t_buf_total * &
                         (T_edges(j+1) - T_edges(j)))
            dP_dlnT(j) = hist_weight(j) / (t_buf_total * &
                         (log(T_edges(j+1)) - log(T_edges(j))))
         end do
      end if

      deallocate(T_buf, dt_buf)
      if (present(n_samples_used)) n_samples_used = n_used
      if (status /= MC_OK .and. .not. present(status_out)) &
         error stop 'mc_run_engine_buffered: incomplete cooling interval'
   end subroutine mc_run_engine_buffered

end module mc_engine
