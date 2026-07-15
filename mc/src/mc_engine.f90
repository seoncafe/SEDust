module mc_engine
   ! Thread-safe MC engine.  Takes an mc_grain_t (filled by one of the
   ! grain_setup_* routines) and an rng_t and produces a P(T) histogram
   ! plus aggregate (t, T) records.  Each thread should own its own
   ! mc_grain_t and rng_t; there is no shared state in this module.

   use constants,      only: wp, pi
   use mc_grain_type,  only: mc_grain_t
   use mc_rng,         only: rng_t, rng_uniform, rng_exp
   use mc_heatcap,     only: U_of_T, C_of_T
   use radfield,       only: J_Mathis, bbody
   use q_graphite_d16_sphere_mod, only: q_graphite_d16_sphere_abs
   implicit none
   private
   public :: grain_setup_graphite, grain_setup_from_cabs
   public :: grain_get_C, grain_get_QPL, grain_T_from_U
   public :: mc_run_engine
   public :: mc_run_engine_2pass, mc_run_engine_buffered
   public :: mc_traj_enable, mc_traj_disable, mc_traj_get

   real(wp), parameter :: c_cgs       = 2.99792458e10_wp
   real(wp), parameter :: h_cgs       = 6.62606957e-27_wp
   real(wp), parameter :: hc_erg_um   = h_cgs * c_cgs * 1.0e4_wp
   real(wp), parameter :: sigma_sb_cgs= 5.6703744e-5_wp

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
   integer,            save :: traj_n       = 0
   integer,            save :: traj_cap     = 0
   real(wp),           save :: traj_stride  = 0.0_wp
   real(wp),           save :: traj_t_last  = -1.0e30_wp
   ! Set by mc_run_engine before each cool_segment call so that
   ! cool_segment can map its internal tau_seg to a global wall time.
   real(wp),           save :: traj_seg_start_time = 0.0_wp
   real(wp), allocatable, save :: traj_time_buf(:), traj_temp_buf(:)

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
   end subroutine mc_traj_enable


   subroutine mc_traj_disable()
      if (allocated(traj_time_buf)) deallocate(traj_time_buf, traj_temp_buf)
      traj_enabled = .false.
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
      if (traj_n >= traj_cap) return
      if (.not. force) then
         if (t_global < traj_t_last + traj_stride) return
      end if
      traj_n = traj_n + 1
      traj_time_buf(traj_n) = t_global
      traj_temp_buf(traj_n) = T_now
      traj_t_last    = t_global
   end subroutine traj_push

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


   subroutine grain_setup_from_cabs(g, a_um, U_isrf, comp, lam_in, Q_in, lam_c_in)
      ! General setup: caller supplies the wavelength grid and Q_abs(lam).
      ! Used by the SED builder for arbitrary grain compositions
      ! (e.g., astrodust grains, PAH grains with xi-blend).
      type(mc_grain_t), intent(out)          :: g
      real(wp),         intent(in)           :: a_um, U_isrf
      character(len=*), intent(in)           :: comp
      real(wp),         intent(in)           :: lam_in(:), Q_in(:)
      real(wp),         intent(in), optional :: lam_c_in

      integer  :: NL, NCDF, NT, i, j
      real(wp), allocatable :: J_SI(:), integrand_cont(:), integrand_event(:)
      real(wp), allocatable :: cum_event(:), BB(:), lam_node(:)
      real(wp) :: w, Fmax, integ_BQ, integ_B, Tt, frac_c, ie_c, ic_c

      NL = size(lam_in)
      if (size(Q_in) /= NL) then
         write(*,'(a)') 'grain_setup_from_cabs: size(lam_in) /= size(Q_in)'
         stop 1
      end if

      g%a_um = a_um
      g%comp = comp
      if (present(lam_c_in)) g%lam_c = lam_c_in

      g%NLAM = NL
      allocate(g%lam_grid(NL), g%Q_grid(NL), g%u_lam_grid(NL))
      g%lam_grid = lam_in
      g%Q_grid   = Q_in

      ! Radiation field u_lam from J_Mathis (SI -> cgs conversion).
      !   u_lam_cgs[erg/cm3/um] = (4 pi / c_SI) * J_SI[W/m2/m/sr] * 1e-5
      allocate(J_SI(NL))
      call J_Mathis(U_isrf, g%lam_grid, J_SI)
      g%u_lam_grid = (4.0_wp * pi / 2.99792458e8_wp) * J_SI * 1.0e-5_wp
      deallocate(J_SI)

      ! Continuous + stochastic split
      allocate(integrand_cont(NL), integrand_event(NL))
      do i = 1, NL
         integrand_cont(i)  = g%Q_grid(i) * g%u_lam_grid(i)
         integrand_event(i) = pi * (a_um*1.0e-4_wp)**2 * g%Q_grid(i) &
                              * c_cgs * g%u_lam_grid(i) * g%lam_grid(i) / hc_erg_um
      end do
      g%H_cont     = 0.0_wp
      g%rate_event = 0.0_wp
      do i = 1, NL-1
         w = 0.5_wp * (g%lam_grid(i+1) - g%lam_grid(i))
         if (g%lam_grid(i+1) <= g%lam_c) then
            ! interval entirely blueward of (or ending at) lam_c: absorption events
            g%rate_event = g%rate_event &
                 + w * (integrand_event(i) + integrand_event(i+1))
         else if (g%lam_grid(i) >= g%lam_c) then
            ! interval entirely redward of (or starting at) lam_c: continuous heating
            g%H_cont = g%H_cont + w * (integrand_cont(i) + integrand_cont(i+1))
         else
            ! interval straddles lam_c: split at lam_c by linear interpolation so
            ! the blueward part feeds rate_event and the redward part feeds H_cont;
            ! otherwise the cutoff bin is dropped from both terms.
            frac_c = (g%lam_c - g%lam_grid(i)) / (g%lam_grid(i+1) - g%lam_grid(i))
            ie_c   = integrand_event(i) + frac_c*(integrand_event(i+1) - integrand_event(i))
            ic_c   = integrand_cont(i)  + frac_c*(integrand_cont(i+1)  - integrand_cont(i))
            g%rate_event = g%rate_event &
                 + 0.5_wp * (g%lam_c - g%lam_grid(i)) * (integrand_event(i) + ie_c)
            g%H_cont = g%H_cont &
                 + 0.5_wp * (g%lam_grid(i+1) - g%lam_c) * (ic_c + integrand_cont(i+1))
         end if
      end do
      g%H_cont = (c_cgs / 4.0_wp) * g%H_cont

      ! A grain with no absorption blueward of lam_c has no stochastic events;
      ! rng_exp would divide by a zero rate.  Such a grain belongs on the
      ! equilibrium path, so stop rather than propagate Inf/NaN into the engine.
      if (g%rate_event <= 0.0_wp) then
         write(*,'(a,es12.4,a)') &
              'grain_setup_from_cabs: grain radius a =', a_um, &
              ' um has no absorption blueward of lam_c; such grains should use the equilibrium path'
         stop 1
      end if

      ! Inverse-CDF table for photon-wavelength sampling
      NCDF = 1024
      g%NCDF = NCDF
      allocate(g%cdf_F(NCDF), g%cdf_lam(NCDF))
      allocate(cum_event(NL), lam_node(NL))
      lam_node     = g%lam_grid
      cum_event(1) = 0.0_wp
      do i = 2, NL
         if (g%lam_grid(i) <= g%lam_c) then
            ! interval entirely blueward of (or ending at) lam_c
            w = 0.5_wp * (g%lam_grid(i) - g%lam_grid(i-1))
            cum_event(i) = cum_event(i-1) + w * (integrand_event(i-1) + integrand_event(i))
         else if (g%lam_grid(i-1) < g%lam_c) then
            ! interval straddles lam_c: accumulate only the [lam(i-1), lam_c] part,
            ! using the same split integral as the rate_event loop so cum_event(NL)
            ! matches rate_event exactly.  Pin this node at lam_c so the inverse CDF
            ! never returns a sampled wavelength redward of the cutoff.
            frac_c = (g%lam_c - g%lam_grid(i-1)) / (g%lam_grid(i) - g%lam_grid(i-1))
            ie_c   = integrand_event(i-1) + frac_c*(integrand_event(i) - integrand_event(i-1))
            cum_event(i) = cum_event(i-1) &
                 + 0.5_wp * (g%lam_c - g%lam_grid(i-1)) * (integrand_event(i-1) + ie_c)
            lam_node(i)  = g%lam_c
         else
            ! interval entirely redward of lam_c: no event mass; pin at lam_c so
            ! inversion in the flat tail cannot land redward of the cutoff.
            cum_event(i) = cum_event(i-1)
            lam_node(i)  = g%lam_c
         end if
      end do
      Fmax = cum_event(NL)
      if (Fmax <= 0.0_wp) Fmax = 1.0_wp
      do j = 1, NCDF
         g%cdf_F(j) = real(j-1, wp) / real(NCDF-1, wp)
      end do
      do j = 1, NCDF
         call invert_cdf(cum_event, lam_node, NL, g%cdf_F(j)*Fmax, g%cdf_lam(j))
      end do
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
      ! Planck-averaged <Q>_T from the pre-tabulated grid by log-T interp.
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


   function grain_T_from_U(g, U_target) result(T)
      ! Invert U(T) given target U_target [erg].
      type(mc_grain_t), intent(in) :: g
      real(wp),         intent(in) :: U_target
      real(wp) :: T, frac
      integer  :: lo, hi, mid
      if (U_target <= g%U_grid(1)) then
         T = max(g%T_grid(1) * sqrt(max(U_target,1.0e-300_wp) / max(g%U_grid(1),1.0e-300_wp)), 0.5_wp)
         return
      end if
      if (U_target >= g%U_grid(g%NT)) then;  T = g%T_grid(g%NT);  return; end if
      lo = 1; hi = g%NT
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (g%U_grid(mid) <= U_target) then;  lo = mid;  else;  hi = mid;  end if
      end do
      if (g%U_grid(hi) > g%U_grid(lo)) then
         frac = (U_target - g%U_grid(lo)) / (g%U_grid(hi) - g%U_grid(lo))
      else
         frac = 0.0_wp
      end if
      T = exp(log(g%T_grid(lo)) + frac * (log(g%T_grid(hi)) - log(g%T_grid(lo))))
   end function grain_T_from_U


   ! =========================================================================
   ! MC engine: integrate cooling, sample events, accumulate P(T)
   ! =========================================================================

   subroutine mc_run_engine(g, rng, T_init, N_events, t_max, &
                            NHIST, T_edges, dP_dT, dP_dlnT, &
                            t_total, e_abs_tot, e_emit_tot, &
                            time_rec_out, temp_pre_rec_out, temp_rec_out, &
                            dtime_rec_out, n_rec_out)
      ! Inputs:
      !   g          : grain state (filled by grain_setup_*)
      !   rng        : RNG state (thread-local)
      !   T_init     : initial T [K]
      !   N_events   : number of stochastic absorption events to simulate
      !   t_max      : maximum simulated time [s]; 0 -> unlimited
      !   NHIST      : number of histogram bins for P(T)
      !   T_edges(NHIST+1)   : log-spaced edges (filled here from 2 K to 5000 K)
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
      integer,          intent(out), optional :: n_rec_out

      real(wp), parameter :: T_HMIN = 2.0_wp, T_HMAX = 5000.0_wp
      real(wp) :: hist_weight(NHIST)
      real(wp) :: Tcur, time_cur, dt_event, dt_act
      real(wp) :: lam_sample, U_before, U_after, T_after, dE_phot
      real(wp) :: T_cool_end, e_emit_step, lT_min, dlnT_bin
      integer  :: ev, j, n_rec, n_rec_cap
      logical  :: do_record

      lT_min   = log(T_HMIN)
      dlnT_bin = (log(T_HMAX) - lT_min) / real(NHIST, wp)
      do j = 1, NHIST+1
         T_edges(j) = exp(lT_min + (j-1)*dlnT_bin)
      end do
      hist_weight = 0.0_wp

      Tcur     = T_init
      time_cur = 0.0_wp
      e_abs_tot  = 0.0_wp
      e_emit_tot = 0.0_wp

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
         dt_event = rng_exp(rng, g%rate_event)
         if (t_max > 0.0_wp .and. time_cur + dt_event > t_max) then
            dt_act = t_max - time_cur
            traj_seg_start_time = time_cur
            call cool_segment(g, Tcur, dt_act, T_cool_end, e_emit_step, &
                              hist_weight, lT_min, dlnT_bin)
            Tcur = T_cool_end
            time_cur = t_max
            e_emit_tot = e_emit_tot + e_emit_step
            if (do_record .and. n_rec < n_rec_cap) then
               n_rec = n_rec + 1
               time_rec_out(n_rec)      = time_cur
               temp_pre_rec_out(n_rec)  = Tcur    ! cooling-only, no jump beyond t_max
               temp_rec_out(n_rec)      = Tcur
               dtime_rec_out(n_rec)     = dt_act
            end if
            exit
         end if
         traj_seg_start_time = time_cur
         call cool_segment(g, Tcur, dt_event, T_cool_end, e_emit_step, &
                           hist_weight, lT_min, dlnT_bin)
         Tcur = T_cool_end
         time_cur = time_cur + dt_event
         e_emit_tot = e_emit_tot + e_emit_step

         lam_sample = sample_photon_lam(g, rng)
         dE_phot    = hc_erg_um / lam_sample
         U_before   = U_of_T(Tcur, g%a_um, g%comp)
         U_after    = U_before + dE_phot
         T_after    = grain_T_from_U(g, U_after)
         e_abs_tot  = e_abs_tot + dE_phot
         Tcur       = T_after
         ! Force-push the post-jump T to capture the vertical event jump.
         call traj_push(time_cur, Tcur, .true.)

         if (do_record .and. n_rec < n_rec_cap) then
            n_rec = n_rec + 1
            time_rec_out(n_rec)      = time_cur
            temp_pre_rec_out(n_rec)  = T_cool_end   ! T just before photon jump
            temp_rec_out(n_rec)      = Tcur         ! T just after photon jump
            dtime_rec_out(n_rec)     = dt_event
         end if
      end do

      if (do_record .and. present(n_rec_out)) n_rec_out = n_rec
      t_total = time_cur

      dP_dT   = 0.0_wp
      dP_dlnT = 0.0_wp
      if (t_total > 0.0_wp) then
         do j = 1, NHIST
            dP_dT(j)   = hist_weight(j) / (t_total * (T_edges(j+1) - T_edges(j)))
            dP_dlnT(j) = hist_weight(j) / (t_total * dlnT_bin)
         end do
      end if
   end subroutine mc_run_engine


   subroutine cool_segment(g, T_in, dt_total, T_out, e_emit, &
                           hist_weight, lT_min, dlnT_bin)
      type(mc_grain_t), intent(in)    :: g
      real(wp),         intent(in)    :: T_in, dt_total, lT_min, dlnT_bin
      real(wp),         intent(out)   :: T_out, e_emit
      real(wp),         intent(inout) :: hist_weight(:)

      real(wp) :: T0, T1, tau_seg, dts, dTdt0, dTdt1, slope, T_eq_loc, a_cm
      real(wp) :: dt_left, T_probe, half_life, rel_excursion
      real(wp) :: e_emit_step_now
      real(wp), parameter :: T_FLOOR    = 0.5_wp
      real(wp), parameter :: DLNT_MAX   = 0.5_wp
      real(wp), parameter :: STEP_FRAC  = 1.0_wp     ! half-life cap (coarse)
      real(wp), parameter :: DT_FRAC_MAX = 0.05_wp   ! max fractional dT per step
                                                     ! (bounds linearization error
                                                     !  of the super-linear cooling)
      real(wp), parameter :: TINY_RATE  = 1.0e-300_wp
      real(wp), parameter :: CONVERGED  = 1.0e-3_wp
      real(wp), parameter :: LN2        = 0.6931471805599453_wp
      integer,  parameter :: MAX_STEPS  = 4000
      integer  :: nstep

      a_cm    = g%a_um * 1.0e-4_wp
      T0      = T_in
      tau_seg = 0.0_wp
      e_emit  = 0.0_wp

      nstep = 0
      ! Trajectory dump: push the starting state of this segment.
      call traj_push(traj_seg_start_time + tau_seg, T0, .true.)
      do while (tau_seg < dt_total .and. nstep < MAX_STEPS)
         nstep   = nstep + 1
         dt_left = dt_total - tau_seg
         dTdt0   = ode_rhs(g, T0, a_cm)
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
               ! Step cap. The exponential advance below is exact only for a
               ! LINEAR ODE; the true cooling is strongly super-linear
               ! (dT/dt propto -Q(T) T^4), so over a large step the local
               ! linearization lags the true (faster) early cooling and the
               ! grain appears to linger at high T -- over-populating the
               ! hot tail of P(T) and over-emitting in the mid-IR (most
               ! severe for the violently heated, tiniest PAH grains).
               ! Bounding the fractional temperature change per step to
               ! DT_FRAC_MAX keeps the linearization error negligible.
               dts = min(dt_left, STEP_FRAC * half_life, &
                         DT_FRAC_MAX * T0 / max(abs(dTdt0), TINY_RATE))
               T1  = T_eq_loc + (T0 - T_eq_loc) * exp(slope * dts)
            end if
         else if (abs(dTdt0) > 0.0_wp) then
            dts = min(dt_left, DT_FRAC_MAX * T0 / abs(dTdt0))
            T1  = T0 + dts * dTdt0
         else
            dts = dt_left
            T1  = T0
         end if
         if (T1 < T_FLOOR) T1 = T_FLOOR

         ! Sub-sample the exponential trajectory to (a) accumulate the
         ! histogram weight of each bin and (b) estimate the emission integral.
         call accumulate_step(g, T0, T1, T_eq_loc, slope, dts, &
                              hist_weight, lT_min, dlnT_bin, a_cm, &
                              e_emit_step_now)
         e_emit = e_emit + e_emit_step_now

         T0      = T1
         tau_seg = tau_seg + dts
         ! Trajectory dump: stride-controlled intermediate sample.
         call traj_push(traj_seg_start_time + tau_seg, T0, .false.)
      end do
      ! Ensure the endpoint is recorded even if stride filter rejected it.
      call traj_push(traj_seg_start_time + tau_seg, T0, .true.)
      T_out = T0
   end subroutine cool_segment


   subroutine accumulate_step(g, T0, T1, T_eq_loc, slope, dts, &
                              hist_weight, lT_min, dlnT_bin, a_cm, e_emit_out)
      ! Sub-sample the cooling step in time, assigning histogram weight to
      ! each bin and accumulating the emission integral over the step.
      ! Returns the step's emission so the caller can apply the
      ! energy-conservation correction.
      type(mc_grain_t), intent(in)    :: g
      real(wp),         intent(in)    :: T0, T1, T_eq_loc, slope, dts
      real(wp),         intent(inout) :: hist_weight(:)
      real(wp),         intent(in)    :: lT_min, dlnT_bin, a_cm
      real(wp),         intent(out)   :: e_emit_out
      integer,  parameter :: N_SUB = 8
      real(wp) :: dts_sub, time_sub, Tval_sub
      integer  :: i, ibin
      dts_sub = dts / real(N_SUB, wp)
      e_emit_out = 0.0_wp
      do i = 1, N_SUB
         time_sub = (real(i,wp) - 0.5_wp) * dts_sub
         if (abs(slope) > 0.0_wp .and. abs(T0 - T_eq_loc) > 0.0_wp) then
            Tval_sub = T_eq_loc + (T0 - T_eq_loc) * exp(slope * time_sub)
         else
            Tval_sub = T0 + (T1 - T0) * (time_sub / dts)
         end if
         if (Tval_sub < 0.5_wp) Tval_sub = 0.5_wp
         ibin = 1 + int((log(Tval_sub) - lT_min) / dlnT_bin)
         if (ibin >= 1 .and. ibin <= size(hist_weight)) then
            hist_weight(ibin) = hist_weight(ibin) + dts_sub
         end if
         e_emit_out = e_emit_out + 4.0_wp * pi * a_cm**2 * grain_get_QPL(g, Tval_sub) &
                      * sigma_sb_cgs * Tval_sub**4 * dts_sub
      end do
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
      type(mc_grain_t), intent(in)    :: g
      type(rng_t),      intent(inout) :: rng
      real(wp) :: lam, u, frac
      integer  :: lo, hi, mid
      u = rng_uniform(rng)
      if (u <= g%cdf_F(1)) then;  lam = g%cdf_lam(1);  return; end if
      if (u >= g%cdf_F(g%NCDF)) then; lam = g%cdf_lam(g%NCDF); return; end if
      lo = 1; hi = g%NCDF
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (g%cdf_F(mid) <= u) then;  lo = mid;  else;  hi = mid;  end if
      end do
      frac = (u - g%cdf_F(lo)) / (g%cdf_F(hi) - g%cdf_F(lo))
      lam  = g%cdf_lam(lo) + frac * (g%cdf_lam(hi) - g%cdf_lam(lo))
   end function sample_photon_lam


   subroutine invert_cdf(cum, x, n, F, x_at_F)
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: cum(n), x(n), F
      real(wp), intent(out) :: x_at_F
      integer  :: lo, hi, mid
      real(wp) :: frac
      if (F <= cum(1)) then;  x_at_F = x(1);  return; end if
      if (F >= cum(n)) then;  x_at_F = x(n);  return; end if
      lo = 1; hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (cum(mid) <= F) then;  lo = mid;  else;  hi = mid;  end if
      end do
      if (cum(hi) > cum(lo)) then
         frac = (F - cum(lo)) / (cum(hi) - cum(lo))
      else
         frac = 0.0_wp
      end if
      x_at_F = x(lo) + frac * (x(hi) - x(lo))
   end subroutine invert_cdf


   ! ==========================================================================
   ! Shared helpers for adaptive-grid engines.
   ! step_advance + sub_sample_step factor out the trajectory math from
   ! cool_segment so that the binning policy (fixed-grid / probe / replay /
   ! buffer) is decoupled from the integration.
   ! ==========================================================================

   subroutine step_advance(g, T0, dt_left, a_cm, T1, dts, T_eq_loc, slope)
      ! One internal cooling step.  Same logic as the inline block in
      ! cool_segment(), exposed so adaptive-grid engines can reuse it.
      type(mc_grain_t), intent(in)  :: g
      real(wp),         intent(in)  :: T0, dt_left, a_cm
      real(wp),         intent(out) :: T1, dts, T_eq_loc, slope
      real(wp) :: dTdt0, dTdt1, T_probe, half_life, rel_excursion
      real(wp), parameter :: T_FLOOR    = 0.5_wp
      real(wp), parameter :: DLNT_MAX   = 0.5_wp
      real(wp), parameter :: STEP_FRAC  = 1.0_wp
      real(wp), parameter :: CONVERGED  = 1.0e-3_wp
      real(wp), parameter :: LN2        = 0.6931471805599453_wp

      dTdt0   = ode_rhs(g, T0, a_cm)
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
            dts = min(dt_left, STEP_FRAC * half_life)
            T1  = T_eq_loc + (T0 - T_eq_loc) * exp(slope * dts)
         end if
      else if (abs(dTdt0) > 0.0_wp) then
         dts      = min(dt_left, DLNT_MAX * T0 / abs(dTdt0))
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
                              is_log, lT_or_T_lo, dlnT_or_dT)
      ! Bin N_SUB equal-dt sub-samples into hist_weight on the adaptive grid.
      integer,  intent(in)    :: N_SUB, NHIST
      real(wp), intent(in)    :: T_sub(N_SUB), dt_sub, lT_or_T_lo, dlnT_or_dT
      logical,  intent(in)    :: is_log
      real(wp), intent(inout) :: hist_weight(NHIST)
      integer  :: i, ibin
      real(wp) :: x
      do i = 1, N_SUB
         if (is_log) then
            x = log(max(T_sub(i), 1.0e-30_wp))
         else
            x = T_sub(i)
         end if
         ibin = 1 + int((x - lT_or_T_lo) / dlnT_or_dT)
         if (ibin >= 1 .and. ibin <= NHIST) then
            hist_weight(ibin) = hist_weight(ibin) + dt_sub
         end if
      end do
   end subroutine bin_sub_samples


   ! ==========================================================================
   ! Method A: two-pass with deterministic RNG replay.
   ! Pass 1 propagates the trajectory and tracks T_min/T_max of all sub-step
   ! samples after a burn-in.  Pass 2 restores the RNG state and re-runs with
   ! the adaptive grid built from the observed range.  Cost ~ 2 x mc_run_engine.
   ! ==========================================================================

   subroutine mc_run_engine_2pass(g, rng, T_init, N_events, N_pass1, N_burn, t_max, &
                                  NHIST, T_edges, dP_dT, dP_dlnT, &
                                  t_total, e_abs_tot, e_emit_tot, &
                                  is_log_out, T_min_obs_out, T_max_obs_out)
      ! Inputs:
      !   N_pass1 : number of events in Pass 1 (the T_min/T_max probe).
      !             Pass 2 always uses the full N_events for the binned
      !             trajectory.  If N_pass1 <= 0, defaults to N_events/2
      !             (the probe converges much faster than the histogram).
      !             If N_pass1 > N_events, silently clamped to N_events.
      !   N_burn  : number of initial events to exclude from min/max tracking.
      !             If <= 0, defaults to max(50, N_pass1_eff/20).  If
      !             >= N_pass1_eff, silently set to 0.
      ! Outputs (optional):
      !   is_log_out, T_min_obs_out, T_max_obs_out : grid diagnostics
      type(mc_grain_t), intent(in)    :: g
      type(rng_t),      intent(inout) :: rng
      real(wp),         intent(in)    :: T_init, t_max
      integer,          intent(in)    :: N_events, N_pass1, N_burn, NHIST
      real(wp),         intent(out)   :: T_edges(NHIST+1)
      real(wp),         intent(out)   :: dP_dT(NHIST), dP_dlnT(NHIST)
      real(wp),         intent(out)   :: t_total, e_abs_tot, e_emit_tot
      logical,          intent(out), optional :: is_log_out
      real(wp),         intent(out), optional :: T_min_obs_out, T_max_obs_out

      integer,  parameter :: N_SUB     = 8
      real(wp), parameter :: MARGIN    = 0.05_wp
      integer,  parameter :: MAX_STEPS = 300

      integer(kind=8) :: rng_saved
      real(wp) :: hist_weight(NHIST)
      real(wp) :: Tcur, time_cur, dt_event, dt_act, a_cm, tau_seg
      real(wp) :: lam_sample, U_before, U_after, dE_phot
      real(wp) :: T_min_obs, T_max_obs
      real(wp) :: T0_step, T1_step, dts, T_eq_loc, slope
      real(wp) :: T_sub(N_SUB), dt_sub, lT_or_T_lo, dlnT_or_dT
      logical  :: is_log
      integer  :: ev, j, nstep, N_burn_eff, N_pass1_eff

      a_cm = g%a_um * 1.0e-4_wp
      N_pass1_eff = N_pass1
      ! Default = N_events (full probe pass).  Shortening Pass 1 makes the
      ! adaptive grid built from N_pass1 < N_events samples too narrow: the
      ! remaining (N_events - N_pass1) samples in Pass 2 that fall outside
      ! the grid are silently dropped by bin_sub_samples, undercounting the
      ! P(T) tails.  The dropped samples are typically at higher T than the
      ! Pass-1 range; the resulting flux deficit shows up most strongly at
      ! the FIR/sub-mm where B_lam(T) on the Rayleigh-Jeans side is linear
      ! in T.  Empirically, N_pass1 = N/2 produces a ~2-10 percentage-point
      ! extra deficit relative to fixed-grid binning at lambda > 80 um.
      if (N_pass1_eff <= 0)         N_pass1_eff = N_events
      if (N_pass1_eff > N_events)   N_pass1_eff = N_events
      N_burn_eff = N_burn
      if (N_burn_eff <= 0)             N_burn_eff = max(50, N_pass1_eff / 20)
      if (N_burn_eff >= N_pass1_eff)   N_burn_eff = 0

      rng_saved = rng%state

      ! ---- Pass 1: probe (only first N_pass1_eff events) ----
      Tcur     = T_init
      time_cur = 0.0_wp
      T_min_obs = huge(1.0_wp)
      T_max_obs = 0.0_wp
      do ev = 1, N_pass1_eff
         dt_event = rng_exp(rng, g%rate_event)
         if (t_max > 0.0_wp .and. time_cur + dt_event > t_max) then
            dt_act = t_max - time_cur
         else
            dt_act = dt_event
         end if

         T0_step = Tcur
         tau_seg = 0.0_wp
         nstep   = 0
         do while (tau_seg < dt_act .and. nstep < MAX_STEPS)
            nstep = nstep + 1
            call step_advance(g, T0_step, dt_act - tau_seg, a_cm, &
                              T1_step, dts, T_eq_loc, slope)
            call sub_sample_step(T0_step, T1_step, T_eq_loc, slope, dts, &
                                 N_SUB, T_sub, dt_sub)
            if (ev > N_burn_eff) then
               T_min_obs = min(T_min_obs, minval(T_sub))
               T_max_obs = max(T_max_obs, maxval(T_sub))
            end if
            T0_step = T1_step
            tau_seg = tau_seg + dts
         end do
         Tcur     = T0_step
         time_cur = time_cur + tau_seg

         if (t_max > 0.0_wp .and. time_cur >= t_max) exit

         lam_sample = sample_photon_lam(g, rng)
         dE_phot    = hc_erg_um / lam_sample
         U_before   = U_of_T(Tcur, g%a_um, g%comp)
         U_after    = U_before + dE_phot
         Tcur       = grain_T_from_U(g, U_after)
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

      ! ---- Pass 2: replay with binning ----
      rng%state   = rng_saved
      Tcur        = T_init
      time_cur    = 0.0_wp
      e_abs_tot   = 0.0_wp
      e_emit_tot  = 0.0_wp
      hist_weight = 0.0_wp

      do ev = 1, N_events
         dt_event = rng_exp(rng, g%rate_event)
         if (t_max > 0.0_wp .and. time_cur + dt_event > t_max) then
            dt_act = t_max - time_cur
         else
            dt_act = dt_event
         end if

         T0_step = Tcur
         tau_seg = 0.0_wp
         nstep   = 0
         do while (tau_seg < dt_act .and. nstep < MAX_STEPS)
            nstep = nstep + 1
            call step_advance(g, T0_step, dt_act - tau_seg, a_cm, &
                              T1_step, dts, T_eq_loc, slope)
            call sub_sample_step(T0_step, T1_step, T_eq_loc, slope, dts, &
                                 N_SUB, T_sub, dt_sub)
            call bin_sub_samples(T_sub, N_SUB, dt_sub, hist_weight, NHIST, &
                                 is_log, lT_or_T_lo, dlnT_or_dT)
            e_emit_tot = e_emit_tot + e_emit_substep(g, T_sub, N_SUB, dt_sub, a_cm)
            T0_step = T1_step
            tau_seg = tau_seg + dts
         end do
         Tcur     = T0_step
         time_cur = time_cur + tau_seg

         if (t_max > 0.0_wp .and. time_cur >= t_max) exit

         lam_sample = sample_photon_lam(g, rng)
         dE_phot    = hc_erg_um / lam_sample
         U_before   = U_of_T(Tcur, g%a_um, g%comp)
         U_after    = U_before + dE_phot
         e_abs_tot  = e_abs_tot + dE_phot
         Tcur       = grain_T_from_U(g, U_after)
      end do

      t_total = time_cur

      dP_dT   = 0.0_wp
      dP_dlnT = 0.0_wp
      if (t_total > 0.0_wp) then
         do j = 1, NHIST
            dP_dT(j)   = hist_weight(j) / (t_total * (T_edges(j+1) - T_edges(j)))
            dP_dlnT(j) = hist_weight(j) / (t_total * &
                         (log(T_edges(j+1)) - log(T_edges(j))))
         end do
      end if
   end subroutine mc_run_engine_2pass


   ! ==========================================================================
   ! Method B: single-pass with sub-step sample buffer.
   ! Trajectory and emission are computed in one pass; sub-step (T, dt) pairs
   ! after the burn-in are appended to caller-bounded arrays.  T_min/T_max
   ! is tracked online (not capped by the buffer).  After the run the
   ! adaptive grid is built and the buffer is binned.
   ! ==========================================================================

   subroutine mc_run_engine_buffered(g, rng, T_init, N_events, N_burn, t_max, &
                                     NHIST, max_samples, &
                                     T_edges, dP_dT, dP_dlnT, &
                                     t_total, e_abs_tot, e_emit_tot, &
                                     is_log_out, T_min_obs_out, T_max_obs_out, &
                                     n_samples_used, buffer_full_out)
      ! Inputs:
      !   max_samples : capacity of internal sub-step (T,dt) buffer.
      !                 Once full, later sub-steps still update T_min/T_max
      !                 but are not added to the histogram, so the histogram
      !                 is built from the FIRST max_samples post-burn samples.
      type(mc_grain_t), intent(in)    :: g
      type(rng_t),      intent(inout) :: rng
      real(wp),         intent(in)    :: T_init, t_max
      integer,          intent(in)    :: N_events, N_burn, NHIST, max_samples
      real(wp),         intent(out)   :: T_edges(NHIST+1)
      real(wp),         intent(out)   :: dP_dT(NHIST), dP_dlnT(NHIST)
      real(wp),         intent(out)   :: t_total, e_abs_tot, e_emit_tot
      logical,          intent(out), optional :: is_log_out, buffer_full_out
      real(wp),         intent(out), optional :: T_min_obs_out, T_max_obs_out
      integer,          intent(out), optional :: n_samples_used

      integer,  parameter :: N_SUB     = 8
      real(wp), parameter :: MARGIN    = 0.05_wp
      integer,  parameter :: MAX_STEPS = 300

      real(wp) :: hist_weight(NHIST)
      real(wp), allocatable :: T_buf(:), dt_buf(:)
      real(wp) :: Tcur, time_cur, dt_event, dt_act, a_cm, tau_seg
      real(wp) :: lam_sample, U_before, U_after, dE_phot
      real(wp) :: T_min_obs, T_max_obs
      real(wp) :: T0_step, T1_step, dts, T_eq_loc, slope
      real(wp) :: T_sub(N_SUB), dt_sub, lT_or_T_lo, dlnT_or_dT, t_buf_total
      logical  :: is_log, buf_full
      integer  :: ev, j, k, nstep, N_burn_eff, n_used

      a_cm = g%a_um * 1.0e-4_wp
      N_burn_eff = N_burn
      if (N_burn_eff <= 0)       N_burn_eff = max(50, N_events / 20)
      if (N_burn_eff >= N_events) N_burn_eff = 0

      allocate(T_buf(max_samples), dt_buf(max_samples))
      n_used   = 0
      buf_full = .false.

      Tcur       = T_init
      time_cur   = 0.0_wp
      e_abs_tot  = 0.0_wp
      e_emit_tot = 0.0_wp
      T_min_obs  = huge(1.0_wp)
      T_max_obs  = 0.0_wp

      do ev = 1, N_events
         dt_event = rng_exp(rng, g%rate_event)
         if (t_max > 0.0_wp .and. time_cur + dt_event > t_max) then
            dt_act = t_max - time_cur
         else
            dt_act = dt_event
         end if

         T0_step = Tcur
         tau_seg = 0.0_wp
         nstep   = 0
         do while (tau_seg < dt_act .and. nstep < MAX_STEPS)
            nstep = nstep + 1
            call step_advance(g, T0_step, dt_act - tau_seg, a_cm, &
                              T1_step, dts, T_eq_loc, slope)
            call sub_sample_step(T0_step, T1_step, T_eq_loc, slope, dts, &
                                 N_SUB, T_sub, dt_sub)
            if (ev > N_burn_eff) then
               T_min_obs = min(T_min_obs, minval(T_sub))
               T_max_obs = max(T_max_obs, maxval(T_sub))
               do k = 1, N_SUB
                  if (n_used >= max_samples) then
                     buf_full = .true.
                     exit
                  end if
                  n_used         = n_used + 1
                  T_buf(n_used)  = T_sub(k)
                  dt_buf(n_used) = dt_sub
               end do
            end if
            e_emit_tot = e_emit_tot + e_emit_substep(g, T_sub, N_SUB, dt_sub, a_cm)
            T0_step = T1_step
            tau_seg = tau_seg + dts
         end do
         Tcur     = T0_step
         time_cur = time_cur + tau_seg

         if (t_max > 0.0_wp .and. time_cur >= t_max) exit

         lam_sample = sample_photon_lam(g, rng)
         dE_phot    = hc_erg_um / lam_sample
         U_before   = U_of_T(Tcur, g%a_um, g%comp)
         U_after    = U_before + dE_phot
         e_abs_tot  = e_abs_tot + dE_phot
         Tcur       = grain_T_from_U(g, U_after)
      end do

      t_total = time_cur

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
                              is_log, lT_or_T_lo, dlnT_or_dT)
         t_buf_total = t_buf_total + dt_buf(k)
      end do

      dP_dT   = 0.0_wp
      dP_dlnT = 0.0_wp
      if (t_buf_total > 0.0_wp) then
         do j = 1, NHIST
            dP_dT(j)   = hist_weight(j) / (t_buf_total * &
                         (T_edges(j+1) - T_edges(j)))
            dP_dlnT(j) = hist_weight(j) / (t_buf_total * &
                         (log(T_edges(j+1)) - log(T_edges(j))))
         end do
      end if

      deallocate(T_buf, dt_buf)
      if (present(n_samples_used)) n_samples_used = n_used
   end subroutine mc_run_engine_buffered

end module mc_engine
