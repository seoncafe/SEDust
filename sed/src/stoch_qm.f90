module stoch_qm_mod
   ! PARAB_INTERPL is the Aigen Li / B. T. Draine routine transcribed from
   ! BDraine_src/src/parab_interpl.f (2000-11-13 history entry). It retains
   ! its original algorithm for state-density interpolation. Radiative
   ! rates instead use the separate nonnegative_rf_cabs scheme.
   !--------------------------------------------------------------------
   ! Quantum-mechanical stochastic heating solver implementing Draine's
   ! 'dbdis' (thermal-discrete) method
   ! (Draine & Li 2001, Draine & Hensley 2021).
   !
   ! Solves for the steady-state probability distribution P(state) over
   ! discrete enthalpy bins, including:
   !   - Heating transitions (photon absorption): upward rates A_{fi}
   !   - Thermal cooling transitions: downward rates A_{fi} using
   !     Planck emission at the upper bin temperature T_i (thermal-
   !     discrete method, DL01 Section 6.1).
   !   - Intrabin transitions (self-transition correction)
   !   - Steady state by GTH state reduction (direct, cancellation-free)
   !   - Emission from P(state) using the thermal-discrete Planck formula:
   !     sum_J P(J) * B_lam(T_J) * Cabs * 8*pi/lam^4
   !
   ! Three cooling treatments are provided:
   !   'dbdis'  thermal-discrete: all downward transitions f < i retained
   !            (the production reference).
   !   'dbcon'  thermal-continuous: the downward rates are collapsed into
   !            the nearest-neighbour transition, reproducing the
   !            continuous-cooling approximation of the Guhathakurta-Draine
   !            temperature-space recursion.
   !   'stati'  exact-statistical (microcanonical): the downward rate
   !            carries the degeneracy ratio g_f/g_i instead of a Planck
   !            factor, with g(U) counted exactly by the Beyer-Swinehart
   !            algorithm over the vibrational mode spectrum.  Emission is
   !            then the degeneracy-weighted kernel rather than a bin
   !            blackbody.  This treatment needs the full mode spectrum and
   !            a feasible state count, so it is used for PAH grains with
   !            a <= 25 A and reverts to 'dbdis' elsewhere.
   !
   ! The bins use the solver's own vibrational-mode enthalpy (see
   ! T_from_U_bisect), not the model's H(T).  A grain is solved only when
   ! the UMIN/UMAX window converges and the spectrum emits the absorbed
   ! power; otherwise qm_solve_grain says why through `status`.
   !
   ! Units: internally CGS (wavelengths cm, energies erg, cross sections
   ! cm^2, radiation field c*u_lambda erg/cm^3/s). The public API
   ! converts from/to our SI/um convention at the entry/exit points.
   !--------------------------------------------------------------------
   use constants, only: wp, pi, c, hc, h_planck, k_boltz, ev2erg, hc_k
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private

   ! Public API
   public :: qm_solve_grain   ! single-grain QM P(state) solver + emission
   public :: qm_nstate_default
   public :: qm_nisrf_max     ! transition-matrix wavelengths (0 = the full field grid)
   public :: qm_status_name   ! name of a qm_solve_grain status code
   public :: QM_SOLVED, QM_INVALID_INPUT, QM_NO_STEADY_STATE, &
             QM_WINDOW_NOT_CONVERGED, QM_ENERGY_NOT_CONSERVED, QM_ATOM_COUNT_CAP, &
             QM_INVALID_RATE, QM_STATUS_MAX
   public :: qm_method        ! 'dbdis' (thermal-discrete, production),
                              ! 'dbcon' (thermal-continuous / GD
                              ! nearest-neighbour collapse), or
                              ! 'stati' (exact-statistical / microcanonical)
   public :: qm_verbose       ! single-grain stderr diagnostics on/off

   ! Method selector (set once before solving). Default = validated dbdis path.
   character(len=5) :: qm_method = 'dbdis'

   ! Diagnostic-output toggle. Default .true. so the CLI drivers keep the
   ! single-grain stderr lines; dust_emission sets it from the model's
   ! `verbose` field so the library path stays silent by default.
   logical, save :: qm_verbose = .true.

   ! Runtime-settable sizes (Draine's production settings: 500 bins and the
   ! full field grid). qm_nisrf_max = 0 uses the full heating field; N >= 2
   ! below the field count coarsens photon absorption to N points, keeping
   ! its tabulated absorbed power. Cooling always uses the full opacity
   ! grid, shared with the emitted spectrum.
   integer :: qm_nstate_default = 200  ! enthalpy bins
   integer :: qm_nisrf_max      = 0    ! heating-field wavelengths

   ! qm_solve_grain status codes
   integer, parameter :: QM_SOLVED               = 0
   integer, parameter :: QM_INVALID_INPUT        = 1
   integer, parameter :: QM_ATOM_COUNT_CAP       = 2  ! explicit modes needed, natom > NATOM_QM_MAX
   integer, parameter :: QM_NO_STEADY_STATE      = 3
   integer, parameter :: QM_WINDOW_NOT_CONVERGED = 4  ! UMIN/UMAX still moving after MAX_ITER
   integer, parameter :: QM_ENERGY_NOT_CONSERVED = 5  ! |P_em/P_abs - 1| > 5%

   integer, parameter :: QM_INVALID_RATE = 6  ! negative or nonfinite physical rate
   integer, parameter :: QM_STATUS_MAX = QM_INVALID_RATE

   ! Physical constants (cgs: hc, h_planck, c, k_boltz, ev2erg, hc_k) are the
   ! exact SI values of the constants module.

   ! Simpson integration tolerance
   real(wp), parameter :: EPS_SIMP = 1.0d-3
   ! Thermal cooling integrals stop at E1 + N_KT*kT (see simpson_therm_cooling)
   real(wp), parameter :: N_KT = 60.0_wp


contains

   !====================================================================
   ! Nonnegative spectral interpolation for radiative transition rates.
   ! Each interval uses a convex combination of its nonnegative endpoint
   ! samples. Unlike a three-point parabola, this is continuous at every
   ! node and cannot create negative absorption or radiation between nodes.
   ! The original Aigen Li / B. T. Draine PARAB_INTERPL
   ! (BDraine_src/src/parab_interpl.f) is not used for rates.
   !====================================================================
   subroutine nonnegative_rf_cabs(wl_cm, nisrf, isrf_wl, isrf, cabs_arr, &
                             crssct, radfld)
      integer, intent(in) :: nisrf
      real(wp), intent(in) :: wl_cm
      real(wp), intent(in) :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out) :: crssct, radfld
      integer :: left, right, mid
      real(wp) :: fraction

      if (wl_cm < isrf_wl(1)) then
         crssct = cabs_arr(1)
         radfld = 0.0_wp
      else if (wl_cm > isrf_wl(nisrf)) then
         crssct = cabs_arr(nisrf) * (isrf_wl(nisrf) / wl_cm)**2
         radfld = 0.0_wp
      else
         left = 1
         right = nisrf
         do while (right - left > 1)
            mid = (left + right) / 2
            if (wl_cm < isrf_wl(mid)) then
               right = mid
            else
               left = mid
            end if
         end do
         fraction = (wl_cm - isrf_wl(left)) / (isrf_wl(right) - isrf_wl(left))
         crssct = (1.0_wp - fraction) * cabs_arr(left) + fraction * cabs_arr(right)
         radfld = (1.0_wp - fraction) * isrf(left) + fraction * isrf(right)
      end if
   end subroutine nonnegative_rf_cabs


   subroutine parab_interpl(x, y, n, t, z)
      implicit none
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: x(n), y(n), t
      real(wp), intent(out) :: z

      integer  :: i, j, k, l, m
      real(wp) :: s

      if (n <= 0) then
         z = 0.0_wp; return
      else if (n == 1) then
         z = y(1); return
      else if (n == 2) then
         z = (y(1)*(t - x(2)) - y(2)*(t - x(1))) / (x(1) - x(2))
         return
      end if

      ! n >= 3: bracket t
      if (t <= x(2)) then
         k = 1; m = 3
      else if (t >= x(n-1)) then
         k = n - 2; m = n
      else
         k = 1; m = n
         do while (iabs(k - m) /= 1)
            l = (k + m) / 2
            if (t < x(l)) then
               m = l
            else
               k = l
            end if
         end do
         if (abs(t - x(k)) < abs(t - x(m))) then
            k = k - 1
         else
            m = m + 1
         end if
      end if

      z = 0.0_wp
      do i = k, m
         s = 1.0_wp
         do j = k, m
            if (j /= i) s = s * (t - x(j)) / (x(i) - x(j))
         end do
         z = z + s * y(i)
      end do
   end subroutine parab_interpl

   !====================================================================
   ! 2b. Smallest spacing between neighbouring field-grid photon energies
   !     hc/lambda_k on [ea, eb] (eb - ea when no grid interval meets it):
   !     the step a quadrature must reach before a zero sum can be read as a
   !     zero integrand.
   !====================================================================
   pure function photon_grid_spacing(ea, eb, nisrf, isrf_wl) result(de_min)
      real(wp), intent(in) :: ea, eb
      integer,  intent(in) :: nisrf
      real(wp), intent(in) :: isrf_wl(nisrf)
      real(wp) :: de_min, e_hi, e_lo
      integer  :: k
      de_min = eb - ea
      do k = 1, nisrf - 1
         e_hi = hc / isrf_wl(k)
         e_lo = hc / isrf_wl(k + 1)
         if (e_hi >= ea .and. e_lo <= eb) de_min = min(de_min, e_hi - e_lo)
      end do
   end function photon_grid_spacing


   !====================================================================
   ! 3. Photon overlap integral G_{fi} between two enthalpy bins.
   !    For cooling (downward) transitions the convention is
   !    initial = upper bin, final = lower bin.
   !====================================================================
   subroutine calc_glu(wl_cm, uia, uib, ufa, ufb, gfi)
      implicit none
      real(wp), intent(in)  :: wl_cm, uia, uib, ufa, ufb
      real(wp), intent(out) :: gfi

      real(wp) :: e_ph, e1, e2, e3, e4, dui, duf

      e_ph = hc / wl_cm

      ! For cooling: initial bin [UIA,UIB] is UPPER, final [UFA,UFB] is LOWER
      ! Energy released = U_initial - U_final, so photon energies span:
      e1 = uia - ufb
      e2 = min(uia - ufa, uib - ufb)
      e3 = max(uia - ufa, uib - ufb)
      e4 = uib - ufa

      dui = uib - uia
      duf = ufb - ufa

      if (e1 /= e2) then
         if (e_ph < e1) then
            gfi = 0.0_wp
         else if (e_ph >= e1 .and. e_ph < e2) then
            gfi = (e_ph - e1) / duf
         else if (e_ph >= e2 .and. e_ph <= e3) then
            if (e2 == e3) then
               gfi = 1.0_wp
            else
               gfi = min(dui, duf) / duf
            end if
         else if (e_ph > e3 .and. e_ph < e4) then
            gfi = (e4 - e_ph) / duf
         else
            gfi = 0.0_wp
         end if
      else
         if (e_ph < e1 .or. e_ph >= e4) then
            gfi = 0.0_wp
         else
            gfi = 1.0_wp
         end if
      end if
   end subroutine calc_glu


   !====================================================================
   ! 4. Heating integrand: HEATING_FUNCT
   !    Heating integrand, following Draine's method.
   !====================================================================
   function heating_funct(e_ph, nisrf, isrf_wl, isrf, cabs_arr, &
                          e1, e2, e3, e4, dui, duf) result(hf)
      implicit none
      real(wp) :: hf
      real(wp), intent(in) :: e_ph
      integer,  intent(in) :: nisrf
      real(wp), intent(in) :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(in) :: e1, e2, e3, e4, dui, duf

      real(wp) :: wl, crssct, radfld, radfld_e, gul_e

      hf = 0.0_wp
      if (e_ph == 0.0_wp) return

      wl = hc / e_ph
      call nonnegative_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      radfld_e = radfld * wl**2 / (h_planck * c)

      if (e1 /= e2) then
         if (e_ph < e1) then
            gul_e = 0.0_wp
         else if (e_ph >= e1 .and. e_ph < e2) then
            gul_e = (e_ph - e1) / dui
         else if (e_ph >= e2 .and. e_ph <= e3) then
            if (e2 == e3) then
               gul_e = 1.0_wp
            else
               gul_e = min(dui, duf) / dui
            end if
         else if (e_ph > e3 .and. e_ph < e4) then
            gul_e = (e4 - e_ph) / dui
         else
            gul_e = 0.0_wp
         end if
         hf = gul_e * crssct * radfld_e
      else
         if (e_ph < e1 .or. e_ph >= e4) then
            hf = 0.0_wp
         else
            hf = crssct * radfld_e
         end if
      end if
   end function heating_funct


   !====================================================================
   ! 5. Heating integrand for last bin (HEATING_FUNCT_N)
   !    Heating integrand (n-branch), following Draine's method.
   !====================================================================
   function heating_funct_n(e_ph, e1, ec, nisrf, isrf_wl, isrf, cabs_arr) result(hf)
      implicit none
      real(wp) :: hf
      real(wp), intent(in) :: e_ph, e1, ec
      integer,  intent(in) :: nisrf
      real(wp), intent(in) :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)

      real(wp) :: wl, crssct, radfld, cradfld_e

      hf = 0.0_wp
      if (e_ph == 0.0_wp) return

      wl = hc / e_ph
      call nonnegative_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      cradfld_e = radfld * wl**2 / (h_planck * c)

      if (e_ph < ec) then
         hf = crssct * cradfld_e * (e_ph - e1) / (ec - e1)
      else
         hf = crssct * cradfld_e
      end if
   end function heating_funct_n


   !====================================================================
   ! 6. Intrabin heating integrand
   !    Intrabin heating integrand, following Draine's method.
   !====================================================================
   function intrabin_heating_funct(e_ph, dui, nisrf, isrf_wl, isrf, cabs_arr) result(hf)
      implicit none
      real(wp) :: hf
      real(wp), intent(in) :: e_ph, dui
      integer,  intent(in) :: nisrf
      real(wp), intent(in) :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)

      real(wp) :: wl, crssct, radfld, cradfld_e

      hf = 0.0_wp
      if (e_ph == 0.0_wp) return

      wl = hc / e_ph
      call nonnegative_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      cradfld_e = radfld * wl**2 / (h_planck * c)
      hf = (1.0_wp - e_ph / dui) * crssct * cradfld_e
   end function intrabin_heating_funct


   !====================================================================
   ! 7. Exact-statistical cooling integrand.
   !    Spontaneous plus stimulated photon emission, with no thermal
   !    (Planck) weight:
   !       G_lu(E) * E^3 * C_abs(E) * (8*pi / (h * (hc)^2))
   !                                 * (1 + (hc)^3/(8*pi*E^3) * u_E)
   !====================================================================
   function stat_cooling_funct(e_ph, nisrf, isrf_wl, isrf, cabs_arr, &
                               e1, e2, e3, e4, dui, duf) result(cf)
      implicit none
      real(wp) :: cf
      real(wp), intent(in) :: e_ph
      integer,  intent(in) :: nisrf
      real(wp), intent(in) :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(in) :: e1, e2, e3, e4, dui, duf

      real(wp) :: wl, crssct, radfld, radfld_e, glu_e

      cf = 0.0_wp
      if (e_ph == 0.0_wp) return

      wl = hc / e_ph
      call nonnegative_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      radfld_e = radfld * wl**2 / (hc * c)

      if (e1 /= e2) then
         if (e_ph < e1) then
            glu_e = 0.0_wp
         else if (e_ph >= e1 .and. e_ph < e2) then
            glu_e = (e_ph - e1) / duf
         else if (e_ph >= e2 .and. e_ph <= e3) then
            if (e2 == e3) then
               glu_e = 1.0_wp
            else
               glu_e = min(dui, duf) / duf
            end if
         else if (e_ph > e3 .and. e_ph < e4) then
            glu_e = (e4 - e_ph) / duf
         else
            glu_e = 0.0_wp
         end if
         cf = glu_e * (e_ph**3) * crssct * &
              (8.0_wp * pi / (h_planck * hc**2)) * &
              (1.0_wp + hc**3 / (8.0_wp * pi * e_ph**3) * radfld_e)
      else
         if (e_ph < e1 .or. e_ph >= e4) then
            cf = 0.0_wp
         else
            cf = (e_ph**3) * crssct * &
                 (8.0_wp * pi / (h_planck * hc**2)) * &
                 (1.0_wp + hc**3 / (8.0_wp * pi * e_ph**3) * radfld_e)
         end if
      end if
   end function stat_cooling_funct


   !====================================================================
   ! 8. Intrabin exact-statistical cooling integrand
   !====================================================================
   function intrabin_stat_cooling_funct(e_ph, dui, nisrf, isrf_wl, isrf, &
                                        cabs_arr) result(cf)
      implicit none
      real(wp) :: cf
      real(wp), intent(in) :: e_ph, dui
      integer,  intent(in) :: nisrf
      real(wp), intent(in) :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)

      real(wp) :: wl, crssct, radfld, radfld_e

      cf = 0.0_wp
      if (e_ph == 0.0_wp) return

      wl = hc / e_ph
      call nonnegative_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      radfld_e = radfld * wl**2 / (hc * c)

      cf = (8.0_wp * pi / (h_planck * hc**2)) * &
           (1.0_wp - e_ph / dui) * (e_ph**3) * crssct * &
           (1.0_wp + hc**3 / (8.0_wp * pi * e_ph**3) * radfld_e)
   end function intrabin_stat_cooling_funct



   !====================================================================
   ! 9. Adaptive Simpson integrators
   !    Generic adaptive Simpson for the heating/cooling integrands.
   !    Each uses an internal function pointer pattern via select case.
   !====================================================================

   ! --- Simpson integrator for interbin heating rate ---
   subroutine simpson_heating_rate(e1, e2, e3, e4, dui, duf, &
                                   nisrf, isrf_wl, isrf, cabs_arr, &
                                   rate_heating)
      implicit none
      integer,  intent(in)  :: nisrf
      real(wp), intent(in)  :: dui, duf
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(inout) :: e1, e2, e3, e4
      real(wp), intent(out) :: rate_heating

      real(wp) :: elow_cut, eupp_cut, hstep, t1, t2, s1, s2, p, xval
      integer  :: n, k

      if (e1 == e4) then
         rate_heating = 0.0_wp; return
      end if

      elow_cut = hc / isrf_wl(nisrf)
      eupp_cut = hc / isrf_wl(1)

      e1 = max(min(e1, eupp_cut), elow_cut)
      e2 = max(min(e2, eupp_cut), elow_cut)
      e3 = max(min(e3, eupp_cut), elow_cut)
      e4 = max(min(e4, eupp_cut), elow_cut)

      if (e1 == e4) then
         rate_heating = 0.0_wp; return
      end if

      n = 1
      hstep = e4 - e1
      t1 = 0.5_wp * hstep * ( &
           heating_funct(e1, nisrf, isrf_wl, isrf, cabs_arr, e1, e2, e3, e4, dui, duf) + &
           heating_funct(e4, nisrf, isrf_wl, isrf, cabs_arr, e1, e2, e3, e4, dui, duf))
      s1 = t1
      do
         p = 0.0_wp
         do k = 0, n - 1
            xval = e1 + (k + 0.5_wp) * hstep
            p = p + heating_funct(xval, nisrf, isrf_wl, isrf, cabs_arr, &
                                  e1, e2, e3, e4, dui, duf)
         end do
         t2 = (t1 + hstep * p) / 2.0_wp
         s2 = (4.0_wp * t2 - t1) / 3.0_wp
         if (abs(s1) > 0.0_wp .and. abs((s2 - s1) / s1) < EPS_SIMP) exit
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) then
            ! A zero sum is accepted only once the samples resolve the
            ! field grid on this interval; a coarser sampling can step over
            ! a band of nonzero field (e.g. a field confined to 4.78-13.6 eV).
            if (hstep <= photon_grid_spacing(e1, e4, nisrf, isrf_wl)) exit
         end if
         t1 = t2; n = n * 2; hstep = hstep / 2.0_wp; s1 = s2
         if (n > 2**16) exit  ! safety limit
      end do
      rate_heating = s2
   end subroutine simpson_heating_rate


   ! --- Simpson integrator for last-bin heating rate ---
   subroutine simpson_heating_rate_n(e1_in, ec, nisrf, isrf_wl, isrf, &
                                     cabs_arr, rate_heating)
      implicit none
      integer,  intent(in)  :: nisrf
      real(wp), intent(in)  :: e1_in, ec
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out) :: rate_heating

      real(wp) :: elow_cut, eupp_cut, emin, hstep, t1, t2, s1, s2, p, xval
      integer  :: n, k

      elow_cut = hc / isrf_wl(nisrf)
      eupp_cut = hc / isrf_wl(1)

      emin = e1_in
      if (emin < elow_cut) emin = elow_cut
      if (emin >= eupp_cut) then
         rate_heating = 0.0_wp; return
      end if

      n = 1
      hstep = eupp_cut - emin
      t1 = 0.5_wp * hstep * ( &
           heating_funct_n(emin,     e1_in, ec, nisrf, isrf_wl, isrf, cabs_arr) + &
           heating_funct_n(eupp_cut, e1_in, ec, nisrf, isrf_wl, isrf, cabs_arr))
      s1 = t1
      do
         p = 0.0_wp
         do k = 0, n - 1
            xval = emin + (k + 0.5_wp) * hstep
            p = p + heating_funct_n(xval, e1_in, ec, nisrf, isrf_wl, isrf, cabs_arr)
         end do
         t2 = (t1 + hstep * p) / 2.0_wp
         s2 = (4.0_wp * t2 - t1) / 3.0_wp
         if (abs(s1) > 0.0_wp .and. abs((s2 - s1) / s1) < EPS_SIMP) exit
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) then
            ! A zero sum is accepted only once the samples resolve the
            ! field grid on this interval; a coarser sampling can step over
            ! a band of nonzero field (e.g. a field confined to 4.78-13.6 eV).
            if (hstep <= photon_grid_spacing(emin, eupp_cut, nisrf, isrf_wl)) exit
         end if
         t1 = t2; n = n * 2; hstep = hstep / 2.0_wp; s1 = s2
         if (n > 2**16) exit
      end do
      rate_heating = s2
   end subroutine simpson_heating_rate_n


   ! --- Simpson integrator for intrabin heating rate ---
   subroutine simpson_intrabin_heating(dui, nisrf, isrf_wl, isrf, cabs_arr, &
                                       rate_heating)
      implicit none
      integer,  intent(in)  :: nisrf
      real(wp), intent(in)  :: dui
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out) :: rate_heating

      real(wp) :: elow_cut, eupp_cut, ea, eb, hstep, t1, t2, s1, s2, p, xval
      integer  :: n, k

      rate_heating = 0.0_wp
      if (dui <= 0.0_wp) return

      elow_cut = hc / isrf_wl(nisrf)
      eupp_cut = hc / isrf_wl(1)
      if (elow_cut > dui) return

      ea = elow_cut
      eb = dui
      if (dui > eupp_cut) eb = eupp_cut

      n = 1
      hstep = eb - ea
      t1 = 0.5_wp * hstep * ( &
           intrabin_heating_funct(ea, dui, nisrf, isrf_wl, isrf, cabs_arr) + &
           intrabin_heating_funct(eb, dui, nisrf, isrf_wl, isrf, cabs_arr))
      s1 = t1
      do
         p = 0.0_wp
         do k = 0, n - 1
            xval = ea + (k + 0.5_wp) * hstep
            p = p + intrabin_heating_funct(xval, dui, nisrf, isrf_wl, isrf, cabs_arr)
         end do
         t2 = (t1 + hstep * p) / 2.0_wp
         s2 = (4.0_wp * t2 - t1) / 3.0_wp
         if (abs(s1) > 0.0_wp .and. abs((s2 - s1) / s1) < EPS_SIMP) exit
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) then
            ! A zero sum is accepted only once the samples resolve the
            ! field grid on this interval; a coarser sampling can step over
            ! a band of nonzero field (e.g. a field confined to 4.78-13.6 eV).
            if (hstep <= photon_grid_spacing(ea, eb, nisrf, isrf_wl)) exit
         end if
         t1 = t2; n = n * 2; hstep = hstep / 2.0_wp; s1 = s2
         if (n > 2**16) exit
      end do
      rate_heating = s2
   end subroutine simpson_intrabin_heating


   ! --- Simpson integrator for the interbin exact-statistical cooling rate ---
   subroutine simpson_stat_cooling(e1_io, e2_io, e3_io, e4_io, dui, duf, &
                                   nisrf, isrf_wl, isrf, cabs_arr, &
                                   rate_cooling)
      implicit none
      integer,  intent(in)    :: nisrf
      real(wp), intent(inout) :: e1_io, e2_io, e3_io, e4_io
      real(wp), intent(in)    :: dui, duf
      real(wp), intent(in)    :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out)   :: rate_cooling

      real(wp) :: elow_cut, eupp_cut, hstep, t1, t2, s1, s2, p, xval
      integer  :: n, k

      elow_cut = hc / isrf_wl(nisrf)
      eupp_cut = hc / isrf_wl(1)

      e1_io = max(min(e1_io, eupp_cut), elow_cut)
      e2_io = max(min(e2_io, eupp_cut), elow_cut)
      e3_io = max(min(e3_io, eupp_cut), elow_cut)
      e4_io = max(min(e4_io, eupp_cut), elow_cut)

      if (e1_io == e4_io) then
         rate_cooling = 0.0_wp; return
      end if

      n = 1
      hstep = e4_io - e1_io
      t1 = 0.5_wp * hstep * ( &
           stat_cooling_funct(e1_io, nisrf, isrf_wl, isrf, cabs_arr, &
                              e1_io, e2_io, e3_io, e4_io, dui, duf) + &
           stat_cooling_funct(e4_io, nisrf, isrf_wl, isrf, cabs_arr, &
                              e1_io, e2_io, e3_io, e4_io, dui, duf))
      s1 = t1
      do
         p = 0.0_wp
         do k = 0, n - 1
            xval = e1_io + (k + 0.5_wp) * hstep
            p = p + stat_cooling_funct(xval, nisrf, isrf_wl, isrf, cabs_arr, &
                                       e1_io, e2_io, e3_io, e4_io, dui, duf)
         end do
         t2 = (t1 + hstep * p) / 2.0_wp
         s2 = (4.0_wp * t2 - t1) / 3.0_wp
         if (abs(s1) > 0.0_wp .and. abs((s2 - s1) / s1) < EPS_SIMP) exit
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) then
            ! A zero sum is accepted only once the samples resolve the
            ! field grid on this interval; a coarser sampling can step over
            ! a band of nonzero field (e.g. a field confined to 4.78-13.6 eV).
            if (hstep <= photon_grid_spacing(e1_io, e4_io, nisrf, isrf_wl)) exit
         end if
         t1 = t2; n = n * 2; hstep = hstep / 2.0_wp; s1 = s2
         if (n > 2**16) exit
      end do
      rate_cooling = s2
   end subroutine simpson_stat_cooling


   ! --- Simpson integrator for the intrabin exact-statistical cooling rate ---
   subroutine simpson_intrabin_stat_cooling(e1_io, e2_io, dui, &
                                            nisrf, isrf_wl, isrf, cabs_arr, &
                                            rate_cooling)
      implicit none
      integer,  intent(in)    :: nisrf
      real(wp), intent(inout) :: e1_io, e2_io
      real(wp), intent(in)    :: dui
      real(wp), intent(in)    :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out)   :: rate_cooling

      real(wp) :: elow_cut, eupp_cut, hstep, t1, t2, s1, s2, p, xval
      integer  :: n, k

      elow_cut = hc / isrf_wl(nisrf)
      eupp_cut = hc / isrf_wl(1)

      e1_io = max(min(e1_io, eupp_cut), elow_cut)
      e2_io = max(min(e2_io, eupp_cut), elow_cut)

      if (e1_io == e2_io) then
         rate_cooling = 0.0_wp; return
      end if

      n = 1
      hstep = e2_io - e1_io
      t1 = 0.5_wp * hstep * ( &
           intrabin_stat_cooling_funct(e1_io, dui, nisrf, isrf_wl, isrf, cabs_arr) + &
           intrabin_stat_cooling_funct(e2_io, dui, nisrf, isrf_wl, isrf, cabs_arr))
      s1 = t1
      do
         p = 0.0_wp
         do k = 0, n - 1
            xval = e1_io + (k + 0.5_wp) * hstep
            p = p + intrabin_stat_cooling_funct(xval, dui, nisrf, isrf_wl, isrf, cabs_arr)
         end do
         t2 = (t1 + hstep * p) / 2.0_wp
         s2 = (4.0_wp * t2 - t1) / 3.0_wp
         if (abs(s1) > 0.0_wp .and. abs((s2 - s1) / s1) < EPS_SIMP) exit
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) then
            ! A zero sum is accepted only once the samples resolve the
            ! field grid on this interval; a coarser sampling can step over
            ! a band of nonzero field (e.g. a field confined to 4.78-13.6 eV).
            if (hstep <= photon_grid_spacing(e1_io, e2_io, nisrf, isrf_wl)) exit
         end if
         t1 = t2; n = n * 2; hstep = hstep / 2.0_wp; s1 = s2
         if (n > 2**16) exit
      end do
      rate_cooling = s2
   end subroutine simpson_intrabin_stat_cooling



   !====================================================================
   ! 9b. Thermal cooling integrand (THERM_COOLING_FUNCT)
   !     Follows Draine's method.
   !
   !     The photon-emission integrand divided by the Planck thermal
   !     factor (exp(E/kT) - 1). This is the 'dbdis' (thermal-discrete)
   !     method of DL01 Section 6.1: the downward transition rate is set
   !     by the thermal Planck emission at temperature T_i of the upper bin,
   !     from the Planck function at the upper bin temperature.
   !
   !     Returns:
   !       G_lu(E) * E^3 * Cabs(E) * (8*pi / (h*(hc)^2))
   !                                * (1 + (hc)^3/(8*pi*E^3) * u_E)
   !                                / (exp(E/(kB*T)) - 1)
   !====================================================================
   function therm_cooling_funct(e_ph, t_upper, nisrf, isrf_wl, isrf, cabs_arr, &
                                e1, e2, e3, e4, dui, duf) result(cf)
      implicit none
      real(wp) :: cf
      real(wp), intent(in) :: e_ph, t_upper
      integer,  intent(in) :: nisrf
      real(wp), intent(in) :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(in) :: e1, e2, e3, e4, dui, duf

      real(wp) :: wl, crssct, radfld, radfld_e, glu_e, boltz_arg

      cf = 0.0_wp
      if (e_ph == 0.0_wp) return
      if (t_upper <= 0.0_wp) return

      boltz_arg = e_ph / (k_boltz * t_upper)
      if (boltz_arg > 500.0_wp) return   ! exp(-500) ~ 0

      wl = hc / e_ph
      call nonnegative_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      radfld_e = radfld * wl**2 / (hc * c)

      ! Evaluate G_lu(E) -- overlap function for two bins
      ! Note: for therm_cooling, G_lu convention uses DUF in denominator
      ! (matching Draine's thermal cooling integrand)
      if (e1 /= e2) then
         if (e_ph < e1) then
            glu_e = 0.0_wp
         else if (e_ph >= e1 .and. e_ph < e2) then
            glu_e = (e_ph - e1) / duf
         else if (e_ph >= e2 .and. e_ph <= e3) then
            if (e2 == e3) then
               glu_e = 1.0_wp
            else
               glu_e = min(dui, duf) / duf
            end if
         else if (e_ph > e3 .and. e_ph < e4) then
            glu_e = (e4 - e_ph) / duf
         else
            glu_e = 0.0_wp
         end if
         cf = glu_e * (e_ph**3) * crssct * &
              (8.0_wp * pi / (h_planck * hc**2)) * &
              (1.0_wp + hc**3 / (8.0_wp * pi * e_ph**3) * radfld_e) / &
              (exp(boltz_arg) - 1.0_wp)
      else
         if (e_ph < e1 .or. e_ph >= e4) then
            cf = 0.0_wp
         else
            cf = (e_ph**3) * crssct * &
                 (8.0_wp * pi / (h_planck * hc**2)) * &
                 (1.0_wp + hc**3 / (8.0_wp * pi * e_ph**3) * radfld_e) / &
                 (exp(boltz_arg) - 1.0_wp)
         end if
      end if
   end function therm_cooling_funct


   !====================================================================
   ! 9c. Intrabin thermal cooling integrand
   !     Follows Draine's method: the intrabin emission integrand
   !     weighted by 1/(exp(E/kT)-1).
   !====================================================================
   function intrabin_therm_cooling_funct(e_ph, dui, t_upper, nisrf, isrf_wl, &
                                          isrf, cabs_arr) result(cf)
      implicit none
      real(wp) :: cf
      real(wp), intent(in) :: e_ph, dui, t_upper
      integer,  intent(in) :: nisrf
      real(wp), intent(in) :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)

      real(wp) :: wl, crssct, radfld, radfld_e, boltz_arg

      cf = 0.0_wp
      if (e_ph == 0.0_wp) return
      if (t_upper <= 0.0_wp) return

      boltz_arg = e_ph / (k_boltz * t_upper)
      if (boltz_arg > 500.0_wp) return

      wl = hc / e_ph
      call nonnegative_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      radfld_e = radfld * wl**2 / (hc * c)

      cf = (8.0_wp * pi / (h_planck * hc**2)) * &
           (1.0_wp - e_ph / dui) * (e_ph**3) * crssct * &
           (1.0_wp + hc**3 / (8.0_wp * pi * e_ph**3) * radfld_e) / &
           (exp(boltz_arg) - 1.0_wp)
   end function intrabin_therm_cooling_funct


   !====================================================================
   ! 9d. Simpson integrator for interbin thermal cooling rate
   !     Simpson integrator for the interbin thermal cooling rate.
   !====================================================================
   subroutine simpson_therm_cooling(e1_io, e2_io, e3_io, e4_io, dui, duf, &
                                    t_upper, nisrf, isrf_wl, isrf, cabs_arr, &
                                    rate_cooling)
      implicit none
      integer,  intent(in)    :: nisrf
      real(wp), intent(inout) :: e1_io, e2_io, e3_io, e4_io
      real(wp), intent(in)    :: dui, duf, t_upper
      real(wp), intent(in)    :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out)   :: rate_cooling

      real(wp) :: hstep, t1, t2, s1, s2, p, xval, ekt, eup
      integer  :: n, k

      ! Short-circuit if E1 > 20*kT (thermal emission negligible)
      if (e1_io > 20.0_wp * k_boltz * t_upper) then
         rate_cooling = 0.0_wp; return
      end if

      if (e1_io == e4_io) then
         rate_cooling = 0.0_wp; return
      end if

      ! The integrand carries 1/(exp(E/kT) - 1) and peaks a few kT above
      ! E1, while [E1, E4] can be hundreds of kT wide (a bin wider than the
      ! photons that cool it).  The range is therefore cut at E1 + N_KT*kT
      ! and the first grid has a step no wider than kT, so the peak is
      ! always sampled.  Truncation error: beyond E1 + N_KT*kT the integrand
      ! falls as x**m exp(-x) (x = E/kT), with m = 3 from E^3, plus 1 from
      ! the linear overlap ramp, plus the power with which Cabs rises toward
      ! short wavelengths (2 in the far infrared).  For m <= 7 the omitted
      ! fraction, Gamma(m+1, x1+N_KT)/Gamma(m+1, x1), is below 1e-17 at
      ! N_KT = 60 for every x1 = E1/kT in [0, 20] (the short-circuit above).
      ekt  = k_boltz * t_upper
      eup  = min(e4_io, e1_io + N_KT * ekt)
      n    = max(1, 2**ceiling(log(max((eup - e1_io) / ekt, 1.0_wp)) / log(2.0_wp)))
      hstep = (eup - e1_io) / real(n, wp)
      t1 = 0.0_wp
      do k = 0, n
         xval = e1_io + real(k, wp) * hstep
         t1 = t1 + merge(0.5_wp, 1.0_wp, k == 0 .or. k == n) * &
              therm_cooling_funct(xval, t_upper, nisrf, isrf_wl, isrf, cabs_arr, &
                                  e1_io, e2_io, e3_io, e4_io, dui, duf)
      end do
      t1 = t1 * hstep
      s1 = t1
      do
         p = 0.0_wp
         do k = 0, n - 1
            xval = e1_io + (k + 0.5_wp) * hstep
            p = p + therm_cooling_funct(xval, t_upper, nisrf, isrf_wl, isrf, cabs_arr, &
                                         e1_io, e2_io, e3_io, e4_io, dui, duf)
         end do
         t2 = (t1 + hstep * p) / 2.0_wp
         s2 = (4.0_wp * t2 - t1) / 3.0_wp
         ! A zero sum on a grid finer than kT is a zero integrand (e.g. no
         ! absorption cross section at these photon energies).
         if (abs(s1) > 0.0_wp .and. abs((s2 - s1) / s1) < EPS_SIMP) exit
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) exit
         t1 = t2; n = n * 2; hstep = hstep / 2.0_wp; s1 = s2
         if (n > 2**16) exit  ! safety limit
      end do
      rate_cooling = s2
   end subroutine simpson_therm_cooling


   !====================================================================
   ! 9e. Simpson integrator for intrabin thermal cooling rate
   !     Simpson integrator for the intrabin thermal cooling rate.
   !====================================================================
   subroutine simpson_intrabin_therm_cooling(e1_io, e2_io, dui, t_upper, &
                                              nisrf, isrf_wl, isrf, cabs_arr, &
                                              rate_cooling)
      implicit none
      integer,  intent(in)    :: nisrf
      real(wp), intent(inout) :: e1_io, e2_io
      real(wp), intent(in)    :: dui, t_upper
      real(wp), intent(in)    :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out)   :: rate_cooling

      real(wp) :: hstep, t1, t2, s1, s2, p, xval, ekt, eup
      integer  :: n, k

      ! Short-circuit if E1 > 20*kT
      if (e1_io > 20.0_wp * k_boltz * t_upper) then
         rate_cooling = 0.0_wp; return
      end if

      if (e1_io == e2_io) then
         rate_cooling = 0.0_wp; return
      end if

      ! Same peak-resolving range and first step as simpson_therm_cooling
      ! (cut at E1 + N_KT*kT, first step <= kT); the extra factor
      ! (1 - E/DUI) <= 1 only lowers the omitted tail.
      ekt  = k_boltz * t_upper
      eup  = min(e2_io, e1_io + N_KT * ekt)
      n    = max(1, 2**ceiling(log(max((eup - e1_io) / ekt, 1.0_wp)) / log(2.0_wp)))
      hstep = (eup - e1_io) / real(n, wp)
      t1 = 0.0_wp
      do k = 0, n
         xval = e1_io + real(k, wp) * hstep
         t1 = t1 + merge(0.5_wp, 1.0_wp, k == 0 .or. k == n) * &
              intrabin_therm_cooling_funct(xval, dui, t_upper, nisrf, isrf_wl, isrf, cabs_arr)
      end do
      t1 = t1 * hstep
      s1 = t1
      do
         p = 0.0_wp
         do k = 0, n - 1
            xval = e1_io + (k + 0.5_wp) * hstep
            p = p + intrabin_therm_cooling_funct(xval, dui, t_upper, nisrf, isrf_wl, isrf, cabs_arr)
         end do
         t2 = (t1 + hstep * p) / 2.0_wp
         s2 = (4.0_wp * t2 - t1) / 3.0_wp
         if (abs(s1) > 0.0_wp .and. abs((s2 - s1) / s1) < EPS_SIMP) exit
         ! A zero sum on a grid finer than kT is a zero integrand.
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) exit
         t1 = t2; n = n * 2; hstep = hstep / 2.0_wp; s1 = s2
         if (n > 2**16) exit
      end do
      rate_cooling = s2
   end subroutine simpson_intrabin_therm_cooling


   !====================================================================
   ! 10. HEATING_AFI -- heating transition rate from bin I to bin F
   !     Heating transition rate, following Draine's method.
   !====================================================================
   subroutine calc_heating_afi(ibin, ui, uia, uib, fbin, uf, ufa, ufb, &
                               nstate, nisrf, isrf_wl, isrf, cabs_arr, &
                               afi_heating)
      implicit none
      integer,  intent(in)  :: ibin, fbin, nstate, nisrf
      real(wp), intent(in)  :: ui, uia, uib, uf, ufa, ufb
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out) :: afi_heating

      real(wp) :: e1, e2, e3, e4, dui, duf, ec, afi_intrabin

      ! Guard: if uf == ui the rate is undefined
      if (abs(uf - ui) < tiny(1.0_wp)) then
         afi_heating = 0.0_wp; return
      end if

      e1 = ufa - uib
      e2 = min(ufa - uia, ufb - uib)
      e3 = max(ufa - uia, ufb - uib)
      e4 = ufb - uia

      dui = uib - uia
      duf = ufb - ufa

      if (fbin /= nstate) then
         call simpson_heating_rate(e1, e2, e3, e4, dui, duf, &
                                   nisrf, isrf_wl, isrf, cabs_arr, afi_heating)
      else
         ! Last bin: all photons with E > UNA - UIB get assigned here
         ec = ufa - uia
         call simpson_heating_rate_n(e1, ec, nisrf, isrf_wl, isrf, &
                                     cabs_arr, afi_heating)
      end if

      ! Intrabin correction for adjacent bins (F = I+1, I /= 1)
      if (fbin == ibin + 1 .and. ibin /= 1) then
         call simpson_intrabin_heating(dui, nisrf, isrf_wl, isrf, cabs_arr, &
                                       afi_intrabin)
         afi_heating = afi_heating + afi_intrabin
      end if

      ! Convert from energy rate to transition rate (s^-1)
      afi_heating = afi_heating / (uf - ui)
   end subroutine calc_heating_afi


   !====================================================================
   ! 11. STAT_COOLING_AFI -- exact-statistical cooling transition rate
   !     from bin I (upper) to bin F (lower).  The microcanonical rate
   !     carries the degeneracy ratio g_f/g_i = exp(lnG_f - lnG_i) in
   !     place of the Planck factor of the thermal treatment.
   !====================================================================
   subroutine calc_stat_cooling_afi(ibin, ui, uia, uib, lngui, &
                                    fbin, uf, ufa, ufb, lnguf, &
                                    nisrf, isrf_wl, isrf, cabs_arr, &
                                    afi_cooling)
      implicit none
      integer,  intent(in)  :: ibin, fbin, nisrf
      real(wp), intent(in)  :: ui, uia, uib, lngui
      real(wp), intent(in)  :: uf, ufa, ufb, lnguf
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out) :: afi_cooling

      real(wp) :: e1, e2, e3, e4, dui, duf
      real(wp) :: afi_cooling0, afi_intrabin_cooling
      real(wp) :: lng_diff, e1c, e2c

      ! Guard: if ui == uf the rate is undefined
      if (abs(ui - uf) < tiny(1.0_wp)) then
         afi_cooling = 0.0_wp; return
      end if

      ! Photon energies for this transition
      e1 = uia - ufb
      e2 = min(uia - ufa, uib - ufb)
      e3 = max(uia - ufa, uib - ufb)
      e4 = uib - ufa

      dui = uib - uia
      duf = ufb - ufa

      lng_diff = lnguf - lngui

      if (ibin /= fbin + 1) then
         ! Non-adjacent downward transition
         call simpson_stat_cooling(e1, e2, e3, e4, dui, duf, &
                                   nisrf, isrf_wl, isrf, cabs_arr, afi_cooling)

         if (lng_diff < -1.0d300 .or. lng_diff > 1.0d300) then
            afi_cooling = 0.0_wp
         else
            afi_cooling = afi_cooling / (ui - uf) * exp(lng_diff)
         end if

      else
         ! Adjacent downward transition (I = F + 1): interbin + intrabin
         call simpson_stat_cooling(e1, e2, e3, e4, dui, duf, &
                                   nisrf, isrf_wl, isrf, cabs_arr, afi_cooling0)

         if (lng_diff < -1.0d300 .or. lng_diff > 1.0d300) then
            afi_cooling0 = 0.0_wp
         else
            afi_cooling0 = afi_cooling0 / (ui - uf) * exp(lng_diff)
         end if

         ! Intrabin cooling correction
         e1c = uia - ufb
         e2c = dui
         call simpson_intrabin_stat_cooling(e1c, e2c, dui, &
                                            nisrf, isrf_wl, isrf, cabs_arr, &
                                            afi_intrabin_cooling)
         afi_cooling = afi_cooling0 + afi_intrabin_cooling / (ui - uf)
      end if
   end subroutine calc_stat_cooling_afi


   !====================================================================
   ! 11b. THERM_COOLING_AFI -- thermal cooling transition rate
   !      from bin I (upper) to bin F (lower), using Planck emission
   !      at temperature T_I.
   !      Thermal cooling transition rate (METHOD='dbdis'), following
   !      Draine's method.
   !
   !      This is the 'dbdis' (thermal-discrete) method of DL01:
   !      the downward rate is determined by the Planck-weighted
   !      spontaneous + stimulated emission at the upper bin's
   !      temperature. This makes the
   !      transition matrix consistent with the thermal-discrete
   !      emission formula already used in qm_emission.
   !====================================================================
   subroutine calc_therm_cooling_afi(ibin, ti, ui, uia, uib, &
                                     fbin, uf, ufa, ufb, &
                                     delta_u2, &
                                     nisrf, isrf_wl, isrf, cabs_arr, &
                                     afi_cooling)
      implicit none
      integer,  intent(in)  :: ibin, fbin, nisrf
      real(wp), intent(in)  :: ti, ui, uia, uib
      real(wp), intent(in)  :: uf, ufa, ufb, delta_u2
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out) :: afi_cooling

      real(wp) :: e1, e2, e3, e4, dui, duf
      real(wp) :: afi_cooling0, afi_intrabin_cooling
      real(wp) :: e1c, e2c

      ! Guard: if ui == uf the rate is undefined
      if (abs(ui - uf) < tiny(1.0_wp)) then
         afi_cooling = 0.0_wp; return
      end if
      if (ti <= 0.0_wp) then
         afi_cooling = 0.0_wp; return
      end if

      ! Photon energies for this transition (same as stat_cooling)
      e1 = uia - ufb
      e2 = min(uia - ufa, uib - ufb)
      e3 = max(uia - ufa, uib - ufb)
      e4 = uib - ufa

      dui = uib - uia
      duf = ufb - ufa

      if (ibin /= fbin + 1) then
         ! Non-adjacent downward transition
         call simpson_therm_cooling(e1, e2, e3, e4, dui, duf, ti, &
                                    nisrf, isrf_wl, isrf, cabs_arr, afi_cooling)

         ! Draine's method: divide by (UI-UF), multiply by DUF/DUI
         ! (or DELTA_U2/DUI for F=1)
         if (fbin /= 1) then
            afi_cooling = afi_cooling / (ui - uf) * duf / dui
         else
            afi_cooling = afi_cooling / (ui - uf) * delta_u2 / dui
         end if

      else
         ! Adjacent downward transition (I = F + 1): interbin + intrabin
         call simpson_therm_cooling(e1, e2, e3, e4, dui, duf, ti, &
                                    nisrf, isrf_wl, isrf, cabs_arr, afi_cooling0)

         if (fbin /= 1) then
            afi_cooling0 = afi_cooling0 / (ui - uf) * duf / dui
         else
            afi_cooling0 = afi_cooling0 / (ui - uf) * delta_u2 / dui
         end if

         ! Intrabin cooling correction (only for F /= 1)
         if (fbin /= 1) then
            e1c = uia - ufb
            e2c = dui
            call simpson_intrabin_therm_cooling(e1c, e2c, dui, ti, &
                                                nisrf, isrf_wl, isrf, cabs_arr, &
                                                afi_intrabin_cooling)
            afi_cooling = afi_cooling0 + afi_intrabin_cooling / (ui - uf)
         else
            afi_cooling = afi_cooling0
         end if
      end if
   end subroutine calc_therm_cooling_afi


   !====================================================================
   ! 12. Build transition matrix (AMATRIX)
   !     Transition-matrix construction, following Draine's method. The
   !     upward rates are shared; `method` selects the downward kernel.
   !====================================================================
   subroutine build_transition_matrix(nstate, u, ua, ub, lngu, t_bins, &
                                      nisrf, isrf_wl, isrf, cabs_arr, &
                                      method, ncool, cooling_wl, cooling_field, cooling_cabs, &
                                      amatrix, negative_count, rates_ok)
      implicit none
      integer,  intent(in)  :: nstate, nisrf, ncool
      real(wp), intent(in) :: cooling_wl(ncool), cooling_field(ncool), cooling_cabs(ncool)
      integer, intent(out) :: negative_count
      logical, intent(out) :: rates_ok
      real(wp), intent(in)  :: u(nstate), ua(nstate), ub(nstate), lngu(nstate)
      real(wp), intent(in)  :: t_bins(nstate)
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      character(len=*), intent(in) :: method   ! 'dbdis', 'dbcon' or 'stati'
                                               ! (thread-private)
      real(wp), intent(out) :: amatrix(nstate, nstate)

      integer  :: ibin, fbin
      real(wp) :: afi_val, adiag, elow_lim, eupp_lim
      real(wp) :: delta_u2

      amatrix = 0.0_wp
      negative_count = 0
      rates_ok = .true.

      ! Pre-compute photon energy limits from the ISRF wavelength range
      elow_lim = hc / isrf_wl(nisrf)   ! lowest photon energy
      eupp_lim = hc / isrf_wl(1)       ! highest photon energy

      ! ---- UPWARD (heating) transitions: ibin -> fbin, fbin > ibin ----
      ! These are identical for all methods (dbdis/dbcon/stati).
      do ibin = 1, nstate
         do fbin = ibin + 1, nstate
            ! Skip if initial bin has zero width (ground state)
            if (ub(ibin) <= ua(ibin) .and. ibin == 1) then
               ! For ground state: treat as delta at U=0
            end if
            ! Minimum photon energy needed: UFA - UIB
            ! Maximum photon energy that contributes: UFB - UIA
            if (ua(fbin) - ub(ibin) > eupp_lim) cycle
            if (ub(fbin) - ua(ibin) < elow_lim .and. fbin /= nstate) cycle
            ! Skip if final bin has zero width
            if (ub(fbin) <= ua(fbin)) cycle

            call calc_heating_afi(ibin, u(ibin), ua(ibin), ub(ibin), &
                                 fbin, u(fbin), ua(fbin), ub(fbin), &
                                 nstate, nisrf, isrf_wl, isrf, cabs_arr, &
                                 afi_val)
            if (afi_val < 0.0_wp) negative_count = negative_count + 1
            amatrix(fbin, ibin) = afi_val
         end do
      end do

      ! ---- DOWNWARD (cooling) transitions: thermal-discrete (dbdis) ----
      if (nstate >= 2) then
         delta_u2 = ub(2) - ua(2)
      else
         delta_u2 = 0.0_wp
      end if

      do ibin = 2, nstate
         if (ub(ibin) <= ua(ibin)) cycle
         ! The T>0 guard only matters for the thermal (dbdis/dbcon) kernel;
         ! the exact-statistical rate carries no temperature.
         if (method /= 'stati' .and. t_bins(ibin) <= 0.0_wp) cycle

         do fbin = ibin - 1, 1, -1
            if (fbin /= 1 .and. ub(fbin) <= ua(fbin)) cycle
            if (u(ibin) <= u(fbin)) cycle

            if (method == 'stati') then
               ! Exact-statistical cooling: the rate follows g_f/g_i and has
               ! no Planck suppression, so the 20kT cutoff must NOT apply.
               call calc_stat_cooling_afi(ibin, u(ibin), ua(ibin), ub(ibin), &
                                          lngu(ibin), &
                                          fbin, u(fbin), ua(fbin), ub(fbin), &
                                          lngu(fbin), &
                                          ncool, cooling_wl, cooling_field, cooling_cabs, &
                                          afi_val)
            else
               ! Thermal-discrete (dbdis): Planck-suppressed beyond ~20kT,
               ! so distant downward transitions are negligible and skipped.
               if (ua(ibin) - ub(fbin) > 20.0_wp * k_boltz * t_bins(ibin)) exit

               call calc_therm_cooling_afi(ibin, t_bins(ibin), &
                                           u(ibin), ua(ibin), ub(ibin), &
                                           fbin, u(fbin), ua(fbin), ub(fbin), &
                                           delta_u2, &
                                           ncool, cooling_wl, cooling_field, cooling_cabs, &
                                           afi_val)
            end if
            if (afi_val < 0.0_wp) negative_count = negative_count + 1
            amatrix(fbin, ibin) = afi_val
         end do
      end do

      ! Inspect every physical rate before continuous-cooling collapse can
      ! hide a negative component. Preserve the count and refuse the solve.
      rates_ok = negative_count == 0 .and. all(ieee_is_finite(amatrix))
      if (.not. rates_ok) return

      ! ---- dbcon collapse: lump all downward rates into nearest neighbour ----
      ! Following Draine's method. For each upper
      ! bin I, the total cooling POWER sum_F A(F,I)*(U_I-U_F) is preserved but
      ! redistributed entirely into the I -> I-1 transition; all longer
      ! downward jumps are zeroed. This is the continuous-cooling (GD-family)
      ! approximation: the mean cooling rate is exact, the dispersion is lost.
      if (method == 'dbcon') then
         do ibin = 2, nstate
            adiag = 0.0_wp                       ! reuse adiag as the power sum
            do fbin = ibin - 1, 1, -1
               adiag = adiag + amatrix(fbin, ibin) * (u(ibin) - u(fbin))
               if (ibin /= fbin + 1) amatrix(fbin, ibin) = 0.0_wp
            end do
            if (u(ibin) > u(ibin - 1)) &
               amatrix(ibin - 1, ibin) = adiag / (u(ibin) - u(ibin - 1))
         end do
      end if

      ! Diagonal: A_{ii} = -sum_{f/=i} A_{fi}
      do ibin = 1, nstate
         adiag = 0.0_wp
         do fbin = 1, nstate
            adiag = adiag + amatrix(fbin, ibin)
         end do
         amatrix(ibin, ibin) = -adiag
      end do
   end subroutine build_transition_matrix


   !====================================================================
   ! 13. Steady state of the transition matrix by GTH state reduction.
   !
   !     AMATRIX(F,I) (F /= I) is the rate of the transition I -> F; the
   !     steady state solves sum_I A(F,I) P(I) = 0 with sum P = 1.  The
   !     Grassmann-Taksar-Heyman elimination (Grassmann, Taksar & Heyman
   !     1985, Operations Research 33, 1107; bibliographic details not
   !     checked against the published article) removes the states from
   !     the top down, redistributing the flow through each removed state
   !     over its exits.  It uses only the off-diagonal rates and adds
   !     nonnegative numbers, so no cancellation occurs however widely the
   !     rates differ (they span ~16 decades in a hard field); Gaussian
   !     elimination on A, or an iterative solver, loses that accuracy.
   !     Cost: n^3/3 multiply-adds and one n x n work array.
   !
   !     Returns ok = .false. when a state that is entered has no exit to
   !     the states below it (no steady state exists on this bin set).
   !====================================================================
   subroutine stationary_distribution_gth(n, amatrix, p, ok)
      implicit none
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: amatrix(n, n)
      real(wp), intent(out) :: p(n)
      logical,  intent(out) :: ok

      real(wp), allocatable :: r(:,:)     ! r(i,j) = rate i -> j
      real(wp) :: s
      real(wp), parameter :: BIG = 1.0e150_wp
      integer  :: i, j, k

      ok = .true.
      allocate(r(n, n))
      do j = 1, n
         do i = 1, n
            if (i == j) then
               r(i, j) = 0.0_wp
            else
               ! A negative physical rate is not assumed to be round-off.
               if (.not. ieee_is_finite(amatrix(j, i)) .or. amatrix(j, i) < 0.0_wp) then
                  ok = .false.
                  p = 0.0_wp
                  deallocate(r)
                  return
               end if
               r(i, j) = amatrix(j, i)
            end if
         end do
         r(j, j) = 0.0_wp
      end do

      do k = n, 2, -1
         s = sum(r(k, 1:k-1))
         if (s <= 0.0_wp) then
            if (any(r(1:k-1, k) > 0.0_wp)) then
               ok = .false.
               p = 0.0_wp
               deallocate(r)
               return
            end if
            cycle            ! never entered: P(k) = 0
         end if
         r(1:k-1, k) = r(1:k-1, k) / s
         do j = 1, k - 1
            if (r(k, j) == 0.0_wp) cycle
            do i = 1, k - 1
               r(i, j) = r(i, j) + r(i, k) * r(k, j)
            end do
         end do
      end do

      ! Back substitution.  P can span more than the double-precision range
      ! across the bins (a window far below the grain's typical enthalpy),
      ! so the running values are rescaled whenever one exceeds BIG; the
      ! states it pushes below the underflow limit carry no probability.
      p(1) = 1.0_wp
      do k = 2, n
         p(k) = sum(p(1:k-1) * r(1:k-1, k))
         if (p(k) > BIG) p(1:k) = p(1:k) / p(k)
      end do
      p = p / sum(p)
      deallocate(r)
   end subroutine stationary_distribution_gth


   !====================================================================
   ! 14. Debye-limit degeneracy ln(g) for an enthalpy bin.
   !     g(U) dU ~ U^(3N-8) dU for 3N-7 harmonic modes in the Debye
   !     (high-U) limit, so ln g_i = (3N-7) ln(Ub_i / Ua_i).
   !
   !     This limit is smooth in U but not monotonic in the bin index once
   !     the bins are mode-tracking, so it is only a labeling for the
   !     thermal kernels (which never read lngu).  The exact-statistical
   !     kernel needs the counted degeneracies of beyer_swinehart instead;
   !     build_enthalpy_bins_qm reports which of the two it produced.
   !====================================================================
   subroutine compute_debye_lngu(nstate, ua, ub, natom, lngu)
      implicit none
      integer,  intent(in)  :: nstate, natom
      real(wp), intent(in)  :: ua(nstate), ub(nstate)
      real(wp), intent(out) :: lngu(nstate)

      integer  :: i
      real(wp) :: nmodes

      nmodes = real(3 * natom - 7, wp)
      lngu(1) = 0.0_wp
      do i = 2, nstate
         if (ua(i) > 0.0_wp .and. ub(i) > ua(i)) then
            lngu(i) = nmodes * log(ub(i) / ua(i))
         else
            lngu(i) = 0.0_wp
         end if
      end do
   end subroutine compute_debye_lngu


   !====================================================================
   ! 14a. Heapsort -- sort array RA(1:N) in ascending order.
   !      Heapsort (Numerical Recipes).
   !====================================================================
   subroutine hpsort(n, ra)
      implicit none
      integer,  intent(in)    :: n
      real(wp), intent(inout) :: ra(n)
      integer  :: i, ir, j, l
      real(wp) :: rra

      if (n < 2) return
      l  = n / 2 + 1
      ir = n
      do
         if (l > 1) then
            l   = l - 1
            rra = ra(l)
         else
            rra    = ra(ir)
            ra(ir) = ra(1)
            ir     = ir - 1
            if (ir == 1) then
               ra(1) = rra
               return
            end if
         end if
         i = l
         j = l + l
         do while (j <= ir)
            if (j < ir) then
               if (ra(j) < ra(j + 1)) j = j + 1
            end if
            if (rra < ra(j)) then
               ra(i) = ra(j)
               i = j
               j = j + j
            else
               j = ir + 1
            end if
         end do
         ra(i) = rra
      end do
   end subroutine hpsort


   !====================================================================
   ! 14b. SIL_VIBRATIONAL_MODES -- vibrational mode spectrum for silicate
   !      Silicate vibrational-mode spectrum, following Draine's method.
   !
   !      Synthetic mode spectrum:
   !        2*(NAT-2) modes: 2-D Debye, Theta = 500 K
   !        (NAT-2)   modes: 3-D Debye, Theta = 1500 K
   !      Modes are sorted in ascending frequency (cm^-1).
   !
   !      Returns:
   !        nmodes        = 3*(natom - 2) vibrational degrees of freedom
   !        emodes(1:nmodes) = mode frequencies in cm^-1 (sorted ascending)
   !====================================================================
   subroutine sil_vibrational_modes(natom, nmodes, emodes)
      implicit none
      integer,  intent(in)  :: natom
      integer,  intent(out) :: nmodes
      real(wp), allocatable, intent(out) :: emodes(:)

      integer  :: j, nm_2d, nm_3d, ntot
      real(wp) :: beta, djoff, em2d, em3d
      real(wp), parameter :: TD2D = 500.0_wp, TD3D = 1500.0_wp

      ! Convert Debye temperatures to cm^-1 (T / (hc/k), hc/k in cm K)
      em2d = TD2D / hc_k
      em3d = TD3D / hc_k

      nm_2d = 2 * (natom - 2)       ! 2/3 of modes
      nm_3d = natom - 2             ! 1/3 of modes
      ntot  = nm_2d + nm_3d         ! = 3*(natom-2)
      nmodes = ntot

      allocate(emodes(ntot))

      ! 2-D Debye modes (out-of-plane analogue for silicate)
      djoff = 0.5_wp
      beta  = (real(nm_2d, wp)**(1.0_wp/3.0_wp) - 1.0_wp) / &
              (2.0_wp * nm_2d - 1.0_wp)

      emodes(1) = sqrt((1.0_wp - beta) * (1.0_wp - djoff) / real(nm_2d, wp) &
                       + beta) * em2d
      emodes(2) = sqrt((1.0_wp - beta) * (2.0_wp - djoff - 0.5_wp) / &
                       real(nm_2d, wp) + beta) * em2d
      emodes(3) = sqrt((1.0_wp - beta) * (3.0_wp - djoff - 0.5_wp) / &
                       real(nm_2d, wp) + beta) * em2d
      do j = 4, nm_2d
         emodes(j) = em2d * sqrt((1.0_wp - beta) * (real(j, wp) - djoff) / &
                                  real(nm_2d, wp) + beta)
      end do

      ! 3-D Debye modes
      beta = 0.0_wp
      emodes(nm_2d + 1) = em3d * ((1.0_wp - beta) * (1.0_wp - djoff) / &
                                   real(nm_3d, wp) + beta)**(1.0_wp/3.0_wp)
      emodes(nm_2d + 2) = em3d * ((1.0_wp - beta) * (2.0_wp - djoff - 0.5_wp) / &
                                   real(nm_3d, wp) + beta)**(1.0_wp/3.0_wp)
      emodes(nm_2d + 3) = em3d * ((1.0_wp - beta) * (3.0_wp - djoff - 0.5_wp) / &
                                   real(nm_3d, wp) + beta)**(1.0_wp/3.0_wp)
      do j = 4, nm_3d
         emodes(nm_2d + j) = em3d * ((1.0_wp - beta) * (real(j, wp) - djoff) / &
                                      real(nm_3d, wp) + beta)**(1.0_wp/3.0_wp)
      end do

      ! Sort modes in ascending frequency
      call hpsort(nmodes, emodes)
   end subroutine sil_vibrational_modes


   !====================================================================
   ! 14c. PAH_SIZE_ATOMS -- compute NC, NH for a PAH grain of radius a_cm
   !      PAH size-to-atom-count relation, following Draine's method.
   !      NC = 460e21 * a^3;  H/C depends on NC.
   !====================================================================
   subroutine pah_size_atoms(a_cm, nc_out, nh_out)
      implicit none
      real(wp), intent(in)  :: a_cm
      integer,  intent(out) :: nc_out, nh_out
      real(wp) :: nc_r, nh_r

      nc_r = 460.0d21 * a_cm**3
      if (nc_r <= 25.0_wp) then
         nh_r = 0.5_wp * nc_r
      else if (nc_r <= 100.0_wp) then
         nh_r = 0.5_wp * sqrt(25.0_wp * nc_r)
      else
         nh_r = 0.25_wp * nc_r
      end if

      nc_out = max(nint(nc_r), 3)
      nh_out = max(nint(nh_r), 0)
   end subroutine pah_size_atoms


   !====================================================================
   ! 14d. PAH_VIBRATIONAL_MODES -- mode spectrum for graphite/PAH grains
   !      PAH vibrational-mode spectrum, following Draine's method.
   !
   !      For PAH with NC carbon atoms and NH hydrogen atoms:
   !        (NC-2)    C-C out-of-plane modes (2-D Debye, 600 cm^-1)
   !        2*(NC-2)  C-C in-plane modes (2-D Debye, 1740 cm^-1)
   !        NH C-H out-of-plane (886 cm^-1)
   !        NH C-H in-plane     (1161 cm^-1)
   !        NH C-H stretching   (3030 cm^-1)
   !      Sorted ascending.
   !
   !      Returns:
   !        nmodes        = 3*(NC + NH) - 6 vibrational modes
   !        emodes(1:nmodes) = mode frequencies in cm^-1 (sorted ascending)
   !====================================================================
   subroutine pah_vibrational_modes(nc, nh, nmodes, emodes)
      implicit none
      integer,  intent(in)  :: nc, nh
      integer,  intent(out) :: nmodes
      real(wp), allocatable, intent(out) :: emodes(:)

      integer  :: j, n1, n_md, n_s, ntot, offset
      real(wp) :: beta, djoff
      real(wp), parameter :: EMCC_IP = 1740.0_wp   ! C-C in-plane Debye freq (cm^-1)
      real(wp), parameter :: EMCC_OP =  600.0_wp   ! C-C out-of-plane Debye freq
      real(wp), parameter :: EMCH_IP = 1161.0_wp   ! C-H in-plane
      real(wp), parameter :: EMCH_OP =  886.0_wp   ! C-H out-of-plane
      real(wp), parameter :: EMCH_ST = 3030.0_wp   ! C-H stretching

      ! Total modes: 3*(NC-2) C-C modes + 3*NH C-H modes
      ! = 3*NC - 6 + 3*NH = 3*(NC+NH) - 6
      ntot = 3 * (nc + nh) - 6
      if (ntot < 1) ntot = 1
      allocate(emodes(ntot))
      emodes = 0.0_wp
      nmodes = 0
      djoff = 0.5_wp

      ! C-C out-of-plane modes: (NC-2) modes, 2-D Debye with EMCC_OP
      n1   = nc - 2
      n_md = 52
      n_s  = 102
      if (nc <= 54) then
         beta = 0.0_wp
      else if (nc > 54 .and. nc < 102) then
         beta = real(n1 - n_md, wp) / real(n_md * (2 * n1 - 1), wp)
      else
         beta = ((real(n1, wp) / n_md) * (real(n_s, wp) / nc)**(2.0_wp/3.0_wp) &
                 - 1.0_wp) / real(2 * n1 - 1, wp)
      end if

      if (n1 >= 1) then
         emodes(1) = sqrt((1.0_wp - beta) * (1.0_wp - djoff) / real(n1, wp) &
                          + beta) * EMCC_OP
      end if
      if (n1 >= 2) then
         emodes(2) = sqrt((1.0_wp - beta) * (1.5_wp - djoff) / real(n1, wp) &
                          + beta) * EMCC_OP
      end if
      if (n1 >= 3) then
         emodes(3) = sqrt((1.0_wp - beta) * (2.5_wp - djoff) / real(n1, wp) &
                          + beta) * EMCC_OP
      end if
      do j = 4, n1
         emodes(j) = sqrt((1.0_wp - beta) * (real(j, wp) - djoff) / &
                          real(n1, wp) + beta) * EMCC_OP
      end do
      nmodes = nmodes + max(n1, 0)

      ! C-C in-plane modes: 2*(NC-2) modes, 2-D Debye with EMCC_IP
      n1 = 2 * (nc - 2)
      offset = nmodes
      if (n1 >= 1) then
         emodes(offset + 1) = sqrt((1.0_wp - beta) * (1.0_wp - djoff) / &
                                    real(n1, wp) + beta) * EMCC_IP
      end if
      if (n1 >= 2) then
         emodes(offset + 2) = sqrt((1.0_wp - beta) * (1.5_wp - djoff) / &
                                    real(n1, wp) + beta) * EMCC_IP
      end if
      if (n1 >= 3) then
         emodes(offset + 3) = sqrt((1.0_wp - beta) * (2.5_wp - djoff) / &
                                    real(n1, wp) + beta) * EMCC_IP
      end if
      do j = 4, n1
         emodes(offset + j) = sqrt((1.0_wp - beta) * (real(j, wp) - djoff) / &
                                    real(n1, wp) + beta) * EMCC_IP
      end do
      nmodes = nmodes + max(n1, 0)

      ! C-H modes (NH of each type)
      if (nh > 0) then
         offset = nmodes
         do j = 1, nh
            emodes(offset + j) = EMCH_OP
         end do
         nmodes = nmodes + nh

         offset = nmodes
         do j = 1, nh
            emodes(offset + j) = EMCH_IP
         end do
         nmodes = nmodes + nh

         offset = nmodes
         do j = 1, nh
            emodes(offset + j) = EMCH_ST
         end do
         nmodes = nmodes + nh
      end if

      ! Sort modes in ascending frequency
      call hpsort(nmodes, emodes)
   end subroutine pah_vibrational_modes


   !====================================================================
   ! 14e. BEYER_SWINEHART -- density of states by exact state counting.
   !
   !      Given NMODES harmonic mode frequencies (cm^-1), counts the number
   !      of vibrational quantum states g(E) on the grid E = j*DE by the
   !      Beyer-Swinehart recursion: adding one mode of quantum dje shifts
   !      and accumulates the running count, so after all modes g(j) is the
   !      exact number of ways to distribute j quanta.
   !
   !      The count grows roughly as exp(sqrt(E)) and overflows IEEE double
   !      well below the top of the enthalpy range, so the grid is truncated
   !      where the count is still representable and continued above that
   !      with the thermodynamic asymptotic d ln g = dE/kT + d ln(dE), the
   !      temperature being the exact mode-sum inverse of E(T).
   !
   !      Returns:
   !        jemax_out = index of the last energy point
   !        energy(0:jemax_out) = energy grid (cm^-1)
   !        lngsum(0:jemax_out) = ln[cumulative number of states up to E]
   !        ebs_out = energy (erg) below which the direct count was feasible
   !====================================================================
   subroutine beyer_swinehart(nmodes, emodes, de, emax_erg, &
                              jemax_out, energy, lngsum, ebs_out)
      implicit none
      integer,  intent(in)    :: nmodes
      real(wp), intent(in)    :: emodes(nmodes), de
      real(wp), intent(in)    :: emax_erg
      integer,  intent(out)   :: jemax_out
      real(wp), allocatable, intent(out) :: energy(:), lngsum(:)
      real(wp), intent(out)   :: ebs_out

      ! Number of grid points reserved for the asymptotic continuation
      integer,  parameter :: NEXT_MAX = 2000
      real(wp), allocatable :: gcount(:), gsum_tmp(:)
      real(wp), allocatable :: lng_ext(:)     ! ln(g) on the same grid
      real(wp) :: emax_cm1, fac, deltae, e_mid, temp, de2
      integer  :: jemax, jm, dje, je, j, ntotal, jebs

      emax_cm1 = emax_erg / hc

      jemax = nint(emax_cm1 / de)
      if (jemax < 1) jemax = 1

      ! Working arrays: counted part [0:jemax] plus room for the continuation
      ntotal = jemax + NEXT_MAX
      allocate(gcount(0:ntotal))
      allocate(gsum_tmp(0:ntotal))
      allocate(energy(0:ntotal))
      allocate(lngsum(0:ntotal))
      allocate(lng_ext(0:ntotal))

      ! Initialize: ground state g(0) = 1, all others 0
100   continue
      gcount(0:jemax) = 0.0_wp
      gcount(0) = 1.0_wp

      ! Beyer-Swinehart recursion over the mode spectrum
      do jm = 1, nmodes
         dje = nint(emodes(jm) / de)
         if (dje == 0) dje = 1   ! do not allow rounding to zero
         do je = dje, jemax
            gcount(je) = gcount(je) + gcount(je - dje)
         end do
         ! Truncate the grid and restart before the count overflows
         if (gcount(jemax) > 1.0d305) then
            if (jm < nmodes / 32) then
               jemax = int(jemax / 2.0_wp)
            else if (jm < nmodes / 16) then
               jemax = int(jemax / 1.4_wp)
            else if (jm < nmodes / 8) then
               jemax = int(jemax / 1.2_wp)
            else if (jm < nmodes / 4) then
               jemax = int(jemax / 1.1_wp)
            else
               jemax = int(jemax / 1.05_wp)
            end if
            if (jemax < 1) jemax = 1
            goto 100
         end if
      end do

      ! Cumulative sum and logarithms over the counted part
      ebs_out = jemax * de * hc

      gsum_tmp(0) = gcount(0)
      do je = 1, jemax
         gsum_tmp(je) = gsum_tmp(je - 1) + gcount(je)
      end do

      do je = 0, jemax
         energy(je)  = real(je, wp) * de
         lng_ext(je) = log(gcount(je) + 1.0d-99)
         lngsum(je)  = log(gsum_tmp(je))
      end do

      ! Continue above the counted range with the dS = dE/T asymptotic.
      ! Three pieces matter and none may be dropped: the d ln(dE) bin-width
      ! term, the widening grid (double DE until the fractional width reaches
      ! 0.005, then geometric with f = 1.01), and the exact mode-sum
      ! temperature.  Without them ln(g) grows incorrectly and drops at the
      ! counted/asymptotic boundary, making the degeneracies non-monotonic
      ! and the exact-statistical rates overflow.
      if (jemax * de < emax_cm1 .and. energy(jemax) > 0.0_wp) then
         jebs = jemax          ! last exactly counted bin

         ! Step 1: double the bin width until the fractional width >= 0.005
         de2 = de
         je  = jemax
         do
            if (de2 / energy(je) >= 0.005_wp) exit
            if (je + 1 > ntotal) exit
            je  = je + 1
            de2 = 2.0_wp * de2
            energy(je) = energy(je - 1) + de2
            if (energy(je) >= emax_cm1) exit
         end do

         ! Step 2: geometric spacing with f = 1.01 up to emax
         if (energy(je) < emax_cm1 .and. je < ntotal) then
            fac = 1.01_wp
            j   = int(log(emax_cm1 / energy(je)) / log(fac)) + 1
            if (j > ntotal - je) j = ntotal - je
            do jm = 1, j
               energy(je + jm) = energy(je) * fac**jm
            end do
            jemax = je + j
         else
            jemax = je
         end if

         ! Fill ln(g), ln(gsum) over the continuation (je = jebs+1 .. jemax)
         do je = jebs + 1, jemax
            e_mid  = 0.5_wp * (energy(je) + energy(je - 1))   ! arithmetic midpoint
            deltae = energy(je) - energy(je - 1)
            ! Exact vibrational mode sum, not the Debye model, so that the
            ! temperature matches the spectrum the count was made from.
            temp = T_from_U_modes(e_mid * hc, nmodes, emodes)
            if (temp <= 0.0_wp) temp = 1.0_wp
            ! d ln g = dE/kT + d ln(dE), with dE in cm^-1 and T in K, so the
            ! coefficient is hc/k = 1.43878 cm K.  Draine's dens_states.f
            ! carries 1.48377 here (a digit transposition; every other use of
            ! hc/k in his source is 1.43877).  This branch is never reached
            ! for the grains the exact-statistical path runs on, so the
            ! correction does not change any output of the released solver.
            lng_ext(je) = lng_ext(je - 1) &
                          + hc_k * deltae / temp &
                          + log(deltae / (energy(je - 1) - energy(je - 2)))
            if ((lng_ext(je) - lngsum(je - 1)) < 100.0_wp) then
               lngsum(je) = lngsum(je - 1) + &
                            log(1.0_wp + exp(lng_ext(je) - lngsum(je - 1)))
            else
               lngsum(je) = lng_ext(je)
            end if
         end do
      end if

      jemax_out = jemax

      deallocate(gcount, gsum_tmp, lng_ext)
   end subroutine beyer_swinehart


   !====================================================================
   ! 14f. T_FROM_U_MODES -- invert E(T) = sum_j E_j/(exp(E_j/kT)-1) for the
   !      temperature of a grain holding energy U_target (erg), using the
   !      explicit mode spectrum (emodes in cm^-1) rather than a Debye
   !      model.  Bisection; returns T in K, and 0 K for U <= 0.
   !====================================================================
   function T_from_U_modes(U_target, nmodes, emodes) result(T)
      implicit none
      real(wp),         intent(in) :: U_target
      integer,          intent(in) :: nmodes
      real(wp),         intent(in) :: emodes(nmodes)
      real(wp) :: T

      real(wp) :: U_cm1, T_lo, T_hi, T_mid, E_mid
      integer  :: iter

      if (U_target <= 0.0_wp .or. nmodes <= 0) then
         T = 0.0_wp; return
      end if

      ! Convert erg -> cm^-1:  1 cm^-1 = hc = 1.98645e-16 erg
      U_cm1 = U_target / hc

      ! Bracket: E(T) is monotonically increasing
      T_lo = 0.5_wp
      T_hi = 5000.0_wp
      do iter = 1, 20
         E_mid = thermal_energy_modes(T_hi, nmodes, emodes)
         if (E_mid >= U_cm1) exit
         T_hi = T_hi * 2.0_wp
      end do

      ! Bisection
      do iter = 1, 60
         T_mid = 0.5_wp * (T_lo + T_hi)
         E_mid = thermal_energy_modes(T_mid, nmodes, emodes)
         if (E_mid < U_cm1) then
            T_lo = T_mid
         else
            T_hi = T_mid
         end if
         if (abs(T_hi - T_lo) < 0.01_wp) exit
      end do
      T = 0.5_wp * (T_lo + T_hi)
   end function T_from_U_modes


   !====================================================================
   ! 14f2. THERMAL_ENERGY_MODES -- E(T) = sum_j E_j / (exp(E_j/kT) - 1)
   !       over the explicit harmonic mode spectrum.  Energies in cm^-1;
   !       hc/k = 1.43878 cm K converts a mode energy in cm^-1 to kT.
   !====================================================================
   function thermal_energy_modes(T, nmodes, emodes) result(Etot)
      implicit none
      real(wp), intent(in) :: T
      integer,  intent(in) :: nmodes
      real(wp), intent(in) :: emodes(nmodes)
      real(wp) :: Etot

      real(wp), parameter :: hc_over_k = hc_k   ! cm K
      real(wp) :: x
      integer  :: j

      Etot = 0.0_wp
      if (T <= 0.0_wp) return
      do j = 1, nmodes
         x = hc_over_k * emodes(j) / T
         if (x < 500.0_wp) then
            Etot = Etot + emodes(j) / (exp(x) - 1.0_wp)
         end if
      end do
   end function thermal_energy_modes


   !====================================================================
   ! 14g. T_FROM_U_BISECT -- invert U(T) -> T for enthalpy bin labeling.
   !      Uses our Debye enthalpy_DL01 via bisection.
   !      U_target is in erg; returned T is in K.
   ! The QM bin temperatures come from the QM solver's own enthalpy
   ! (vibrational-mode counting, enthalpy_DL01), not from the model's
   ! continuous H(T) used by the other solvers. The resulting difference from
   ! the heuristic/draine solvers (e.g. astrodust S2) is intentional and is NOT
   ! a solver defect. Do not change this enthalpy to the model H(T).
   !====================================================================
   function T_from_U_bisect(U_target, grain_type, a_cm) result(T)
      use enthalpy, only: enthalpy_DL01
      implicit none
      real(wp),         intent(in) :: U_target, a_cm
      character(len=*), intent(in) :: grain_type
      real(wp) :: T

      real(wp) :: T_lo, T_hi, T_mid, U_mid, R_um
      integer  :: iter
      character(len=4) :: dtype

      if (U_target <= 0.0_wp) then
         T = 0.0_wp; return
      end if

      R_um = a_cm * 1.0d4   ! cm -> um
      if (trim(grain_type) == 'sil') then
         dtype = 'Sil '
      else
         dtype = 'Car0'
      end if

      T_lo = 1.0_wp
      T_hi = 5000.0_wp
      do iter = 1, 20
         U_mid = enthalpy_DL01(T_hi, R_um, dtype)
         if (U_mid >= U_target) exit
         T_hi = T_hi * 2.0_wp
      end do

      do iter = 1, 60
         T_mid = 0.5_wp * (T_lo + T_hi)
         U_mid = enthalpy_DL01(T_mid, R_um, dtype)
         if (U_mid < U_target) then
            T_lo = T_mid
         else
            T_hi = T_mid
         end if
         if (abs(T_hi - T_lo) < 0.01_wp) exit
      end do
      T = 0.5_wp * (T_lo + T_hi)
   end function T_from_U_bisect


   !====================================================================
   ! 14h. BUILD_ENTHALPY_BINS_QM -- mode-aware enthalpy bin construction
   !      Enthalpy-bin construction, following Draine's method.
   !
   !      Key feature: when UMIN = 0, bins 2 through NSET track individual
   !      vibrational mode energies (the critical fine structure at low
   !      enthalpy that makes the transition matrix well-conditioned).
   !      Bins NSET+1..NSTATE transition from linear to log spacing.
   !
   !      Arguments:
   !        grain_type = 'sil' or 'pah'
   !        a_cm       = grain radius (cm)
   !        natom      = number of atoms
   !        nc, nh     = carbon/hydrogen atom counts (PAH only; 0 for sil)
   !        umin       = min enthalpy (erg), 0 => stochastic from ground state
   !        umax       = max enthalpy (erg)
   !        nstate     = requested number of bins
   !      Returns:
   !        u, ua, ub   = bin center, lower bound, upper bound (erg)
   !        t_bins      = temperature (K) labeling each bin
   !        lngu        = ln(degeneracy) of each bin
   !        nset_out    = number of mode-tracking bins used
   !        lngu_counted = .true. iff lngu came from the exact
   !                       Beyer-Swinehart count, .false. if it is the
   !                       Debye limit (which the 'stati' kernel may not use)
   !====================================================================
   subroutine build_enthalpy_bins_qm(grain_type, a_cm, natom, nc, nh, &
                                      umin, umax, nstate, &
                                      u, ua, ub, t_bins, lngu, nset_out, &
                                      lngu_counted)
      implicit none
      character(len=*), intent(in)  :: grain_type
      real(wp),         intent(in)  :: a_cm, umin, umax
      integer,          intent(in)  :: natom, nc, nh, nstate
      real(wp),         intent(out) :: u(nstate), ua(nstate), ub(nstate)
      real(wp),         intent(out) :: t_bins(nstate), lngu(nstate)
      integer,          intent(out) :: nset_out
      logical,          intent(out) :: lngu_counted

      integer  :: nmodes, i, nset, nset_max, jset
      real(wp) :: de, dlgu, umax_cm1, ebs
      integer  :: jemax_bs, jemax_est
      real(wp), allocatable :: emodes(:), emodes_new(:)
      real(wp), allocatable :: energy_bs(:), lngsum_bs(:)
      logical  :: skip_modes

      ! ---- Step 1: Get vibrational mode spectrum ----
      ! With UMIN /= 0 the bins are purely log-spaced and the exact
      ! degeneracy count is out of reach anyway, so the explicit mode
      ! spectrum is unnecessary and is skipped.  It is needed only for the
      ! mode-tracking bins of the UMIN = 0 branch below.  Skipping keeps very
      ! large grains (natom > 1e6) feasible.
      skip_modes = (umin > 0.0_wp)
      if (skip_modes) then
         nmodes = 0
         allocate(emodes(1));  emodes = 0.0_wp
      else if (trim(grain_type) == 'sil') then
         call sil_vibrational_modes(natom, nmodes, emodes)
      else
         call pah_vibrational_modes(nc, nh, nmodes, emodes)
      end if

      ! ---- Step 2: Determine DE (bin-boundary snapping quantum, and the
      !      energy grid spacing of the Beyer-Swinehart count) ----
      de = 1.0_wp
      if (nmodes >= 2) then
         call select_de(emodes(2) - emodes(1), de)
      end if

      ! ---- Step 3: Build enthalpy bins ----
      if (umin == 0.0_wp) then
         ! Mode-tracking bins: Draine's convention
         ! Bin 1 = ground state (zero width)
         u(1)  = 0.0_wp
         ua(1) = 0.0_wp
         ub(1) = 0.0_wp

         ! Construct EMODES_NEW = distinct mode energies in sorted order
         allocate(emodes_new(min(nmodes, 1000)))
         nset_max = 0
         do i = 1, nmodes - 1
            if (emodes(i + 1) /= emodes(i)) then
               nset_max = nset_max + 1
               emodes_new(nset_max) = emodes(i)
               if (nset_max >= 1000) exit
            end if
         end do

         ! NSET = number of mode-tracking bins (Draine's method caps this at 11)
         nset = (nset_max + 4) / 2
         if (nset >= 11) nset = 11
         if (nset > nstate - 2) nset = max(nstate - 2, 2)
         nset_out = nset

         ! Safety: make sure we have enough distinct modes
         if (nset_max < 3 .or. nset < 3) then
            ! Too few distinct mode energies to track; bin U logarithmically
            nset = 0; nset_out = 0
            call log_spaced_energy_bins(nstate, umin, umax, u, ua, ub)
            deallocate(emodes_new)
            goto 500   ! jump to the T(U) and degeneracy computation
         end if

         ! Bin 2: centered on first distinct mode energy
         ua(2) = (3.0_wp * emodes_new(1) - emodes_new(2)) / 2.0_wp
         ub(2) = (emodes_new(1) + emodes_new(2)) / 2.0_wp
         u(2)  = (ua(2) + ub(2)) / 2.0_wp

         ! Bin 3: centered on second distinct mode energy
         ua(3) = ub(2)
         ub(3) = (emodes_new(2) + emodes_new(3)) / 2.0_wp
         u(3)  = (ua(3) + ub(3)) / 2.0_wp

         ! Bins 4 - NSET: each contain two normal modes
         do i = 4, nset
            if (i * 2 - 3 > nset_max) then
               ! Not enough distinct modes; cap NSET here
               nset = i - 1; nset_out = nset; exit
            end if
            ub(i) = (emodes_new(i * 2 - 4) + emodes_new(i * 2 - 3)) / 2.0_wp
            ua(i) = ub(i - 1)
            u(i)  = (ua(i) + ub(i)) / 2.0_wp
         end do

         ! Beyond NSET: first linear then log-spaced to UMAX
         umax_cm1 = umax / hc
         jset = nstate   ! fallback

         do i = nset + 1, nstate
            u(i) = u(nset) + (ub(nset) - ua(nset)) * real(i - nset, wp)
            if (i == nstate) then
               jset = i
               exit
            end if
            dlgu = log10(umax_cm1 / u(i)) / real(nstate - i, wp)
            if ((u(i) * (10.0_wp**dlgu) - u(i)) > (u(i) - u(i - 1))) then
               jset = i
               exit
            end if
         end do

         ! Fill linear portion
         do i = nset + 1, jset
            u(i) = u(nset) + (ub(nset) - ua(nset)) * real(i - nset, wp)
         end do

         ! Fill log-spaced portion
         if (jset < nstate) then
            dlgu = log10(umax_cm1 / u(jset)) / real(nstate - jset, wp)
            do i = jset + 1, nstate
               u(i) = u(jset) * 10.0_wp**(real(i - jset, wp) * dlgu)
            end do
         end if

         ! Set bin boundaries for upper bins
         do i = nset + 1, nstate - 1
            ua(i) = (u(i) + u(i - 1)) / 2.0_wp
            ub(i) = (u(i) + u(i + 1)) / 2.0_wp
         end do
         ua(nstate) = ub(nstate - 1)
         ub(nstate) = u(nstate)

         ! Snap bin boundaries to multiples of DE (Draine's convention)
         ! and enforce monotonicity
         do i = 1, nstate
            ua(i) = de * real(nint(ua(i) / de), wp)
            ub(i) = de * real(nint(ub(i) / de), wp)
         end do
         do i = 2, nstate
            ua(i) = ub(i - 1)
            if (ub(i) < ua(i) + de) ub(i) = ua(i) + de
         end do
         do i = 2, nstate
            u(i) = 0.5_wp * (ua(i) + ub(i))
         end do

         ! Convert from cm^-1 to erg
         do i = 1, nstate
            u(i)  = u(i)  * hc
            ua(i) = ua(i) * hc
            ub(i) = ub(i) * hc
         end do

         deallocate(emodes_new)

      else
         ! UMIN /= 0: log-space all bins (already in erg)
         nset = 0; nset_out = 0
         dlgu = log10(umax / umin) / real(nstate - 1, wp)
         do i = 1, nstate
            u(i) = umin * 10.0_wp**(real(i - 1, wp) * dlgu)
         end do
         do i = 2, nstate - 1
            ua(i) = (u(i) + u(i - 1)) / 2.0_wp
            ub(i) = (u(i) + u(i + 1)) / 2.0_wp
         end do
         ua(1) = u(1)
         ub(1) = ua(2)
         ua(nstate) = ub(nstate - 1)
         ub(nstate) = u(nstate)
      end if

      ! ---- Step 4: Compute T(U) for each bin ----
      ! The QM bin temperatures come from the QM solver's own enthalpy
      ! (vibrational-mode counting, enthalpy_DL01), not from the model's
      ! continuous H(T) used by the other solvers. The resulting difference from
      ! the heuristic/draine solvers (e.g. astrodust S2) is intentional and is NOT
      ! a solver defect. Do not change this enthalpy to the model H(T).
500   continue
      do i = 1, nstate
         t_bins(i) = T_from_U_bisect(u(i), grain_type, a_cm)
      end do

      ! ---- Step 5: Degeneracy ln(g) of each bin ----
      ! The exact count is O(nmodes * jemax) in time and O(jemax) in memory,
      ! so it is attempted only when both stay bounded:
      !   (a) grid size: UMAX / (HC * DE) < 5e5  (memory)
      !   (b) work: nmodes * jemax < 5e7          (time)
      ! Otherwise the bins keep their mode-tracking structure but carry the
      ! Debye-limit degeneracy, which the 'stati' kernel may not use.
      if (de > 0.0_wp) then
         jemax_est = nint(umax / (hc * de))
      else
         jemax_est = 0
      end if

      if (.not. skip_modes .and. jemax_est > 0 .and. jemax_est < 500000 .and. &
          real(nmodes, wp) * real(jemax_est, wp) < 5.0d7) then
         call beyer_swinehart(nmodes, emodes, de, umax, jemax_bs, &
                              energy_bs, lngsum_bs, ebs)

         lngu(1) = 0.0_wp
         do i = 2, nstate
            call lngu_from_counted_states(ua(i), ub(i), ebs, de, &
                                          jemax_bs, energy_bs, lngsum_bs, &
                                          lngu(i))
         end do
         if (allocated(energy_bs)) deallocate(energy_bs)
         if (allocated(lngsum_bs)) deallocate(lngsum_bs)
         lngu_counted = .true.
      else
         call compute_debye_lngu(nstate, ua, ub, natom, lngu)
         lngu_counted = .false.
      end if

      if (allocated(emodes)) deallocate(emodes)

   contains

      ! Select DE to resolve mode spacing
      subroutine select_de(dx, de_out)
         real(wp), intent(in)  :: dx
         real(wp), intent(out) :: de_out
         if (dx > 2.0_wp) then
            de_out = 1.0_wp
         else if (dx > 1.0_wp) then
            de_out = 0.5_wp
         else if (dx > 0.4_wp) then
            de_out = 0.2_wp
         else if (dx > 0.2_wp) then
            de_out = 0.1_wp
         else if (dx > 0.1_wp) then
            de_out = 0.05_wp
         else if (dx > 0.04_wp) then
            de_out = 0.02_wp
         else if (dx > 0.02_wp) then
            de_out = 0.01_wp
         else if (dx > 0.01_wp) then
            de_out = 5.0d-3
         else if (dx > 0.004_wp) then
            de_out = 2.0d-3
         else if (dx > 0.002_wp) then
            de_out = 1.0d-3
         else if (dx > 0.001_wp) then
            de_out = 5.0d-4
         else if (dx > 0.0004_wp) then
            de_out = 2.0d-4
         else if (dx > 0.0002_wp) then
            de_out = 1.0d-4
         else if (dx > 0.0001_wp) then
            de_out = 5.0d-5
         else
            de_out = 1.0d-5
         end if
      end subroutine select_de


      ! Logarithmically spaced bins in internal energy U, spanning [umin, umax].
      ! Used when the vibrational spectrum has too few distinct mode energies for
      ! the mode-tracking bin construction above to resolve them.
      subroutine log_spaced_energy_bins(ns, umin_f, umax_f, u_f, ua_f, ub_f)
         integer,  intent(in)  :: ns
         real(wp), intent(in)  :: umin_f, umax_f
         real(wp), intent(out) :: u_f(ns), ua_f(ns), ub_f(ns)
         integer  :: ii
         real(wp) :: dlg, umin_eff

         umin_eff = max(umin_f, umax_f * 1.0d-8)
         dlg = log10(umax_f / umin_eff) / real(ns - 1, wp)
         do ii = 1, ns
            u_f(ii) = umin_eff * 10.0_wp**(real(ii - 1, wp) * dlg)
         end do
         do ii = 2, ns - 1
            ua_f(ii) = (u_f(ii) + u_f(ii - 1)) / 2.0_wp
            ub_f(ii) = (u_f(ii) + u_f(ii + 1)) / 2.0_wp
         end do
         ua_f(1) = u_f(1) * 0.5_wp
         ub_f(1) = ua_f(2)
         ua_f(ns) = ub_f(ns - 1)
         ub_f(ns) = u_f(ns)
      end subroutine log_spaced_energy_bins


      ! Number of states in one bin [ua_erg, ub_erg] from the cumulative
      ! count: ln g = ln[gsum(UB) - gsum(UA)].  Below the exactly counted
      ! energy EBS the grid is uniform in DE and the lookup is an index;
      ! above it the grid widens, so interpolate.
      subroutine lngu_from_counted_states(ua_erg, ub_erg, ebs_val, de_val, &
                                          jemax_val, energy_val, lngsum_val, &
                                          lngu_val)
         real(wp), intent(in)  :: ua_erg, ub_erg, ebs_val, de_val
         integer,  intent(in)  :: jemax_val
         real(wp), intent(in)  :: energy_val(0:jemax_val), lngsum_val(0:jemax_val)
         real(wp), intent(out) :: lngu_val

         integer  :: nua_idx, nub_idx
         real(wp) :: lngsumua, lngsumub

         if (ua_erg <= ebs_val) then
            nua_idx = nint(ua_erg / (hc * de_val))
            nua_idx = max(0, min(nua_idx, jemax_val))
            lngsumua = lngsum_val(nua_idx)
         else
            call parab_interpl(energy_val(1:jemax_val+1), &
                              lngsum_val(1:jemax_val+1), &
                              jemax_val + 1, ua_erg / hc, lngsumua)
         end if

         if (ub_erg <= ebs_val) then
            nub_idx = nint(ub_erg / (hc * de_val))
            nub_idx = max(0, min(nub_idx, jemax_val))
            lngsumub = lngsum_val(nub_idx)
         else
            call parab_interpl(energy_val(1:jemax_val+1), &
                              lngsum_val(1:jemax_val+1), &
                              jemax_val + 1, ub_erg / hc, lngsumub)
         end if

         if ((lngsumub - lngsumua) >= 1.0d-10) then
            lngu_val = lngsumub + log(1.0_wp - exp(lngsumua - lngsumub))
         else
            lngu_val = -100.0_wp   ! sentinel: bin holds essentially no state
         end if
      end subroutine lngu_from_counted_states

   end subroutine build_enthalpy_bins_qm


   !====================================================================
   ! 15. Compute NATOM from grain radius and composition
   !     Follows Draine's method.
   !====================================================================
   function compute_natom(a_cm, grain_type) result(natom)
      implicit none
      real(wp),         intent(in) :: a_cm
      character(len=*), intent(in) :: grain_type
      integer :: natom

      real(wp) :: nc_real

      ! 'gra' and 'pah' are ONE material here: the carbon count of the grain,
      ! which is what Draine & Li's carbonaceous heat capacity and mode
      ! spectrum are both written in.  'gra' used to fall to the default and
      ! take the SILICATE atom count -- a graphite population solved with
      ! 'qm' (the Zubko model has one, and so does MRN) was then given the
      ! wrong number of vibrational degrees of freedom.
      select case (trim(grain_type))
      case ('pah', 'gra')
         nc_real = 472.0d21 * a_cm**3
         natom = max(nint(nc_real), 3)
      case default
         natom = nint((4.0_wp * pi / 3.0_wp) * a_cm**3 * &
                      3.5_wp * 7.0_wp / (172.0_wp * 1.66d-24))
      end select
      natom = max(natom, 3)
   end function compute_natom


   !====================================================================
   ! 17. QM emission from P(state).
   !
   !     Thermal ('dbdis' / 'dbcon') kernel:
   !       EMISSION(I) = sum_J P(J) * (8*pi*hcc/lambda^4) * Cabs
   !                     / (exp(hc/(kT_J * lambda)) - 1)
   !     for bins J where h*nu <= UB(J) (photon energy < bin energy).
   !
   !     This is the thermal-discrete approximation of DL01 Section 6.1:
   !     each bin emits as a blackbody at its temperature T(J), which is
   !     the temperature corresponding to the bin's representative energy
   !     via the T(E) relation of DL01 eq. (32).
   !
   !     Exact-statistical ('stati') kernel: no bin temperature and no
   !     Planck factor.  Each downward transition J -> K carries the
   !     degeneracy ratio g_K/g_J and the bin-overlap weight G_lu, plus the
   !     intrabin term for photons softer than the bin width.  This must be
   !     paired with the matching cooling rates in the transition matrix.
   !
   !     NOTE: The (1+J/B0) stimulated-emission factor is omitted from both
   !     kernels (HD23 publishes net emission only).
   !====================================================================
   subroutine qm_emission(nstate, nisrf, isrf_wl_cm, cabs_cm2, &
                           ua, ub, lngu, pstate, t_bins, &
                           method, emission)
      implicit none
      integer,  intent(in)  :: nstate, nisrf
      real(wp), intent(in)  :: isrf_wl_cm(nisrf), cabs_cm2(nisrf)
      real(wp), intent(in)  :: ua(nstate), ub(nstate)
      real(wp), intent(in)  :: lngu(nstate), pstate(nstate)
      real(wp), intent(in)  :: t_bins(nstate)
      character(len=*), intent(in) :: method   ! 'dbdis', 'dbcon' or 'stati'
                                               ! (thread-private)
      real(wp), intent(out) :: emission(nisrf)

      integer  :: i, j, k
      real(wp) :: hnu, term, hcc_val, dj, pref, gkj

      hcc_val = h_planck * c * c   ! h * c^2

      emission = 0.0_wp

      do i = 1, nisrf
         hnu = hc / isrf_wl_cm(i)

         if (method == 'stati') then
            ! Exact-statistical kernel: degeneracy-weighted, no Planck factor.
            do j = 2, nstate
               if (pstate(j) <= 0.0_wp) cycle
               pref = (8.0_wp * pi * hcc_val / isrf_wl_cm(i)**4) &
                      * cabs_cm2(i) * pstate(j)
               ! Intrabin emission (photon leaves the grain within bin j)
               if ((ub(j) - ua(j)) > 0.0_wp .and. hnu <= ub(j) - ua(j)) then
                  emission(i) = emission(i) + pref &
                                * (1.0_wp - hnu / (ub(j) - ua(j)))
               end if
               ! Interbin emission to all lower bins k < j.  Skip the
               ! sentinel degeneracy (lngu = -100 marks a bin with no state)
               ! and guard the exponential against overflow.
               if (lngu(j) <= -99.0_wp) cycle
               do k = j - 1, 1, -1
                  if (lngu(k) <= -99.0_wp) cycle
                  if (lngu(k) - lngu(j) > 100.0_wp) cycle
                  if (hnu > ub(j) - ua(k)) cycle
                  if (hnu < ua(j) - ub(k)) cycle
                  call calc_glu(isrf_wl_cm(i), ua(j), ub(j), ua(k), ub(k), gkj)
                  if (gkj <= 0.0_wp) cycle
                  emission(i) = emission(i) + pref &
                                * exp(lngu(k) - lngu(j)) * gkj
               end do
            end do

         else
            ! Thermal-discrete (dbdis) kernel: blackbody at bin temperature T_j.
            do j = 2, nstate
               if (pstate(j) <= 0.0_wp) cycle
               if (t_bins(j) <= 0.0_wp) cycle
               ! Photon energy must not exceed the upper bound of this bin
               if (hnu > ub(j)) cycle

               term = hc / (isrf_wl_cm(i) * k_boltz * t_bins(j))
               if (term < 500.0_wp) then
                  dj = (8.0_wp * pi * hcc_val / isrf_wl_cm(i)**4) * &
                       cabs_cm2(i) * pstate(j) / (exp(term) - 1.0_wp)
                  emission(i) = emission(i) + dj
               end if
            end do
         end if
      end do
   end subroutine qm_emission


   !====================================================================
   ! 17b. COARSEN_ABSORBED_FIELD -- a coarser wavelength grid for the
   !      transition matrix that keeps the absorbed power.
   !
   !      NC nodes equally spaced in ln(lambda) span the same range as the
   !      fine grid.  Every fine point (trapezoid weight w_i in lambda)
   !      belongs to its nearest node.  A node carries the interval mean
   !      cross section Cbar = sum C w / sum w and the absorption-weighted
   !      field Jbar = sum C J w / sum C w, so Cbar*Jbar*sum(w) equals the
   !      power absorbed in the interval; the field is then scaled by one
   !      factor so that the trapezoid integral of C*J over the coarse grid
   !      equals that over the fine grid.  The total absorbed power is kept
   !      exactly; the photon spectrum inside an interval is replaced by
   !      its mean, so structure narrower than an interval (an absorption
   !      edge, a band limit of the field) is smoothed over it.  A node
   !      that no fine point falls into takes the fine grid interpolated
   !      linearly in ln(lambda).
   !====================================================================
   subroutine coarsen_absorbed_field(n, wl, fld, cab, nc, wl_c, fld_c, cab_c)
      implicit none
      integer,  intent(in)  :: n, nc
      real(wp), intent(in)  :: wl(n), fld(n), cab(n)
      real(wp), intent(out) :: wl_c(nc), fld_c(nc), cab_c(nc)

      real(wp) :: x1, dx, w, sw(nc), scw(nc), scjw(nc), sjw(nc)
      real(wp) :: p_fine, p_coarse, f
      integer  :: i, k

      x1 = log(wl(1))
      dx = (log(wl(n)) - x1) / real(nc - 1, wp)
      do k = 1, nc
         wl_c(k) = exp(x1 + real(k - 1, wp) * dx)
      end do
      wl_c(1) = wl(1);  wl_c(nc) = wl(n)

      sw = 0.0_wp; scw = 0.0_wp; scjw = 0.0_wp; sjw = 0.0_wp
      do i = 1, n
         w = 0.5_wp * (wl(min(i+1, n)) - wl(max(i-1, 1)))
         k = 1 + nint((log(wl(i)) - x1) / dx)
         k = max(1, min(k, nc))
         sw(k)   = sw(k)   + w
         scw(k)  = scw(k)  + cab(i) * w
         scjw(k) = scjw(k) + cab(i) * fld(i) * w
         sjw(k)  = sjw(k)  + fld(i) * w
      end do
      do k = 1, nc
         if (sw(k) > 0.0_wp) then
            cab_c(k) = scw(k) / sw(k)
            if (scw(k) > 0.0_wp) then
               fld_c(k) = scjw(k) / scw(k)
            else
               fld_c(k) = sjw(k) / sw(k)
            end if
         else
            call linear_in_lnlambda(wl_c(k), cab_c(k), fld_c(k))
         end if
      end do

      p_fine = 0.0_wp
      do i = 1, n - 1
         p_fine = p_fine + 0.5_wp * (cab(i)*fld(i) + cab(i+1)*fld(i+1)) * (wl(i+1) - wl(i))
      end do
      p_coarse = 0.0_wp
      do k = 1, nc - 1
         p_coarse = p_coarse + (wl_c(k+1) - wl_c(k)) / 6.0_wp * &
              (2.0_wp*cab_c(k)*fld_c(k) + cab_c(k)*fld_c(k+1) + &
               cab_c(k+1)*fld_c(k) + 2.0_wp*cab_c(k+1)*fld_c(k+1))
      end do
      if (p_coarse > 0.0_wp) then
         f = p_fine / p_coarse
         fld_c = fld_c * f
      end if

   contains

      subroutine linear_in_lnlambda(x, cx, fx)
         real(wp), intent(in)  :: x
         real(wp), intent(out) :: cx, fx
         integer  :: j
         real(wp) :: t
         j = 1
         do while (j < n - 1 .and. wl(j + 1) < x)
            j = j + 1
         end do
         t = (log(x) - log(wl(j))) / (log(wl(j + 1)) - log(wl(j)))
         t = max(0.0_wp, min(1.0_wp, t))
         cx = (1.0_wp - t) * cab(j) + t * cab(j + 1)
         fx = (1.0_wp - t) * fld(j) + t * fld(j + 1)
      end subroutine linear_in_lnlambda
   end subroutine coarsen_absorbed_field


   !====================================================================
   ! 18. Top-level single-grain QM driver
   !
   ! Given a grain with:
   !   - cross section Cabs(lambda) in cm^2
   !   - radiation field J_lam in SI W/m^3/sr at wavelengths lam_um
   !   - equilibrium enthalpy EEQ (erg), from the caller's model H(T)
   !   - grain type ('sil' or 'pah'; 'gra' takes the carbonaceous modes)
   !   - grain radius a_cm
   ! QM takes no enthalpy argument by design; its bin temperatures come only
   ! from its own mode-counting enthalpy (see T_from_U_bisect). Do not add
   ! one. The resulting difference from the heuristic/draine solvers, which
   ! use the model's continuous H(T) (e.g. astrodust S2), is intentional and
   ! is NOT a solver defect.
   !
   ! Returns emission(NLAM) = lambda * P_lambda in erg/s for one grain,
   ! matching Draine's spectral power convention, and solved = .true.;
   ! otherwise solved = .false., a zero spectrum, and the reason in
   ! `status` (QM_* below), which the caller reports.
   !====================================================================
   subroutine qm_solve_grain(nlam, lam_um, cabs_cm2, j_lam_si, &
                              eeq, eeqss, &
                              a_cm, grain_type, emission_out, solved, status, &
                              actual_cooling_out, raw_negative_count_out)
      ! The single-photon bound is a property of the radiation field, so it is
      ! computed by the radiation-field module and shared with the two call
      ! sites in sed_astrodust_mod rather than re-derived here.
      use radfield,    only: hardest_photon_energy
      implicit none
      integer,          intent(in)  :: nlam
      real(wp),         intent(in)  :: lam_um(nlam), cabs_cm2(nlam)
      real(wp),         intent(in)  :: j_lam_si(nlam)
      real(wp),         intent(in)  :: eeq, eeqss
      real(wp),         intent(in)  :: a_cm
      character(len=*), intent(in)  :: grain_type
      real(wp),         intent(out) :: emission_out(nlam)
      logical,          intent(out) :: solved
      integer, optional, intent(out) :: status
      character(len=5), optional, intent(out) :: actual_cooling_out
      integer, optional, intent(out) :: raw_negative_count_out

      integer :: nstate, nstate1, istat, natom, nc_pah, nh_pah
      integer :: i, iter, jcut, jpmax
      integer :: nset_qm, negative_count, negative_total
      real(wp) :: umin, umax, umaxmin, umaxhi, umaxlo, uminhi, uminlo
      real(wp) :: u_photon_max, pmax, sum_p, ratio_em
      logical  :: refine, lngu_counted, gth_ok, converged, rates_ok
      character(len=5) :: eff_method   ! size-gated requested method
      character(len=5) :: cur_method   ! method actually used this iteration

      ! Grains for which the exact state count is feasible.  Draine limits
      ! 'stati' to a <= 25 Angstrom; above that the count overflows before
      ! it reaches the top of the enthalpy range and the degeneracies stop
      ! being usable, so the thermal kernel takes over.
      real(wp), parameter :: A_STATI_MAX_CM = 25.0d-8   ! 25 Angstrom

      ! Full-resolution CGS arrays for emission computation
      real(wp), allocatable :: isrf_wl_full(:), isrf_full(:), cabs_full(:)
      ! Transition-matrix wavelength grid (the full one unless qm_nisrf_max
      ! asks for fewer points)
      integer :: nisrf_ds
      real(wp), allocatable :: isrf_wl_ds(:), isrf_ds(:), cabs_ds(:)

      real(wp), allocatable :: u(:), ua(:), ub(:), t_bins(:), lngu(:)
      real(wp), allocatable :: amatrix(:,:), pstate(:)

      real(wp), parameter :: PMIN_LO_QM = 1.0d-13
      real(wp), parameter :: PMIN_UP_QM = 1.0d-13
      ! Passes of the UMIN/UMAX refinement.  Each pass moves UMIN at most
      ! 80% of the way to the first populated bin, so a window that starts
      ! at U = 0 for a grain whose population sits two decades higher needs
      ! ~14 passes (measured: astrodust qm_dbcon, nstate = 500); 10 left
      ! such grains unconverged.
      integer,  parameter :: MAX_ITER = 30

      ! Largest atom count for which the explicit-mode treatment (PAH_MODES /
      ! sil modes, ~3*natom modes) is feasible. Above this the grain is
      ! returned unsolved (QM_ATOM_COUNT_CAP) so the caller's GD solver
      ! handles it -- such large grains sit near the equilibrium boundary
      ! anyway. Prevents the ~3e7-element allocation seen for the DL07
      ! carbonaceous population, which extends to large a with non-zero dn.
      integer, parameter :: NATOM_QM_MAX = 1000000

      solved = .false.
      negative_total = 0
      if (present(raw_negative_count_out)) raw_negative_count_out = 0
      if (present(actual_cooling_out)) actual_cooling_out = ''
      emission_out = 0.0_wp
      istat = QM_INVALID_INPUT
      if (present(status)) status = istat
      ratio_em = -1.0_wp
      nset_qm = 0
      natom = 0
      iter = 0
      cur_method = qm_method
      ! Mode-aligned bins need three low-energy bins and two upper bins.
      ! Two field samples are the minimum for the index interpolation below.
      ! qm_nisrf_max = 0 selects the full field grid.
      if (qm_nstate_default < 5 .or. qm_nisrf_max == 1 .or. qm_nisrf_max < 0 &
          .or. nlam < 2) then
         if (qm_verbose) write(*,'(a,3(i0,1x))') &
            'stoch_qm: invalid nstate, nisrf, nlam: ', &
            qm_nstate_default, qm_nisrf_max, nlam
         return
      end if

      if (.not. all(ieee_is_finite(lam_um)) .or. .not. all(ieee_is_finite(cabs_cm2)) .or. &
          .not. all(ieee_is_finite(j_lam_si)) .or. .not. ieee_is_finite(a_cm)) return
      if (any(lam_um <= 0.0_wp) .or. any(cabs_cm2 < 0.0_wp) .or. &
          any(j_lam_si < 0.0_wp) .or. a_cm <= 0.0_wp) return
      if (any(lam_um(2:) <= lam_um(:nlam-1))) return

      select case (qm_method)
      case ('dbdis', 'dbcon', 'stati')
         ! recognized
      case default
         write(*,'(a,a,a)') ' stoch_qm: unknown qm_method "', trim(qm_method), &
                            '" -- expected dbdis, dbcon or stati'
         stop 1
      end select

      ! Which cooling treatment this grain can actually carry.  The
      ! exact-statistical kernel needs the counted degeneracies, which are
      ! reachable only for small PAH grains: for silicate the count that
      ! survives the feasibility guards still leaves a residual
      ! non-monotonicity in ln g (~14 in the log) that exp(lngu(k)-lngu(j))
      ! amplifies into a spurious MIR spike, so silicate uses the thermal
      ! kernel -- as Draine's production does, his silicate 'stati' branch
      ! being flagged untested. Agreement between statistical and thermal
      ! spectra is not assumed, particularly at low excitation.
      eff_method = qm_method
      if (qm_method == 'stati') then
         if (trim(grain_type) /= 'pah' .or. a_cm > A_STATI_MAX_CM) &
            eff_method = 'dbdis'
      end if

      nstate  = qm_nstate_default
      nstate1 = nstate - 1

      ! Number of atoms in the grain
      natom = compute_natom(a_cm, grain_type)

      ! NC, NH of the carbonaceous grain, which the mode spectrum builder
      ! needs; pah_vibrational_modes is the carbonaceous spectrum, graphite
      ! included, so 'gra' asks for it on the same terms as 'pah'.
      nc_pah = 0; nh_pah = 0
      if (trim(grain_type) == 'pah' .or. trim(grain_type) == 'gra') then
         call pah_size_atoms(a_cm, nc_pah, nh_pah)
      end if

      ! Feasibility cap: explicit-mode construction (~3*natom modes) is only
      ! needed when UMIN = 0 (mode-aligned first bins; eeq < 0.1*eeqss) or
      ! for the exact-statistical kernel, which counts states over the mode
      ! spectrum.  With UMIN = eeq/5 > 0 the thermal path has log-spaced bins
      ! and skips the modes, so large grains stay feasible there.
      if (max(natom, nc_pah) > NATOM_QM_MAX .and. &
          (eeq < 0.1_wp * eeqss .or. eff_method == 'stati')) then
         istat = QM_ATOM_COUNT_CAP
         if (present(status)) status = istat
         if (qm_verbose) write(0,'(a,es9.2,a,i0,a,a)') &
            '  QM grain: a[um]=', a_cm * 1.0d4, ' natom=', natom, &
            ' ok=F status=', trim(qm_status_name(istat))
         return
      end if

      ! Convert radiation field from SI to the CGS convention used here.
      ! Draine's method uses c * u_lambda (erg cm^{-3} s^{-1}) at wavelengths in cm.
      ! Our J_lam_SI has dimension [W / m^3 / sr] per unit wavelength (per m).
      !   J_lam_CGS(erg/s/cm^2/sr/cm) = J_lam_SI * 10
      ! because 1 W = 1e7 erg/s, 1 m^-2 = 1e-4 cm^-2, and per m -> per cm
      ! contributes another 1e-2. Therefore:
      !   c*u_lambda = 4*pi * J_lam_CGS = 4*pi * 10 * J_lam_SI
      allocate(isrf_wl_full(nlam), isrf_full(nlam), cabs_full(nlam))
      do i = 1, nlam
         isrf_wl_full(i) = lam_um(i) * 1.0d-4     ! um -> cm
         isrf_full(i)    = 4.0_wp * pi * 10.0_wp * j_lam_si(i)
         cabs_full(i)    = cabs_cm2(i)
      end do

      ! Coarsening applies only to photon absorption. Cooling and emission
      ! share the full opacity grid: absorbed-power preservation does not
      ! preserve the Planck-weighted cooling integral across PAH features.
      if (qm_nisrf_max >= 2 .and. qm_nisrf_max < nlam) then
         nisrf_ds = qm_nisrf_max
         allocate(isrf_wl_ds(nisrf_ds), isrf_ds(nisrf_ds), cabs_ds(nisrf_ds))
         call coarsen_absorbed_field(nlam, isrf_wl_full, isrf_full, cabs_full, &
                                     nisrf_ds, isrf_wl_ds, isrf_ds, cabs_ds)
      else
         nisrf_ds = nlam
         allocate(isrf_wl_ds(nisrf_ds), isrf_ds(nisrf_ds), cabs_ds(nisrf_ds))
         isrf_wl_ds = isrf_wl_full
         isrf_ds    = isrf_full
         cabs_ds    = cabs_full
      end if

      ! Initial UMIN/UMAX guesses (following Draine's method).  The
      ! single-photon term must be the hardest photon the radiation FIELD
      ! carries, i.e. hc/lambda at the shortest wavelength where j_lam_si is
      ! significantly nonzero -- not the short end of the caller's wavelength
      ! grid, which is the grid of the model's optics tables.  The two agree
      ! for astrodust and DL07, whose tables stop at the Lyman limit
      ! (hc/0.0912 um = 13.595 eV < 13.6 eV, so the max() leaves the constant
      ! selected); they differ by a factor 91 for Zubko/ZDA, whose ZDA optics
      ! tables start at 1.0e-3 um while the field carries nothing below the
      ! Lyman limit.  The fixed 13.6 eV is right only for a field that stops
      ! at the Lyman limit; a field carried into the EUV raises the bound.
      u_photon_max = hardest_photon_energy(lam_um, j_lam_si)
      umaxmin = 13.65_wp * ev2erg
      umax = max(max(13.6_wp * ev2erg, u_photon_max) + 2.0_wp * eeq, umaxmin)
      if (eeq < 0.1_wp * eeqss) then
         umin = 0.0_wp
      else
         umin = eeq / 5.0_wp
      end if
      umaxhi = 1.0d70 * hc
      umaxlo = 0.0_wp
      uminhi = 1.0d70 * hc
      uminlo = 0.0_wp

      allocate(u(nstate), ua(nstate), ub(nstate), t_bins(nstate), lngu(nstate))
      allocate(amatrix(nstate, nstate), pstate(nstate))
      pstate = 0.0_wp

      ! Iterative UMIN/UMAX refinement (following Draine's method).  Only a
      ! bin set the refinement accepts gives a solution; a window still
      ! moving after MAX_ITER passes is reported, not used.
      converged = .false.
      do iter = 1, MAX_ITER

         ! Build enthalpy bins with mode-aware QM construction
         call build_enthalpy_bins_qm(grain_type, a_cm, natom, nc_pah, nh_pah, &
                                     umin, umax, nstate, &
                                     u, ua, ub, t_bins, lngu, nset_qm, &
                                     lngu_counted)

         ! The exact-statistical kernel needs counted degeneracies; when the
         ! count is out of reach, this grain uses the thermal one.
         cur_method = eff_method
         if (eff_method == 'stati' .and. .not. lngu_counted) cur_method = 'dbdis'

         call build_transition_matrix(nstate, u, ua, ub, lngu, t_bins, &
                                      nisrf_ds, isrf_wl_ds, isrf_ds, cabs_ds, &
                                      cur_method, nlam, isrf_wl_full, isrf_full, cabs_full, &
                                      amatrix, negative_count, rates_ok)
         negative_total = negative_total + negative_count
         if (.not. rates_ok) then
            istat = QM_INVALID_RATE
            exit
         end if

         call stationary_distribution_gth(nstate, amatrix, pstate, gth_ok)
         if (.not. gth_ok) then
            istat = QM_NO_STEADY_STATE
            exit
         end if

         ! The window tests read the bins below the top one, as Draine's do.
         pmax = maxval(pstate(1:nstate1))
         if (pmax <= 0.0_wp) then
            istat = QM_NO_STEADY_STATE
            exit
         end if
         jpmax = maxloc(pstate(1:nstate1), dim=1)

         ! Adjust UMAX
         refine = .false.

         if (pstate(nstate1) / pmax <= PMIN_UP_QM .and. umax > umaxmin) then
            umaxhi = umax
            jcut = nstate1
            do i = nstate1, 1, -1
               jcut = i
               if (ub(i) < umaxmin) exit
               if (pstate(i) / pmax > PMIN_UP_QM) exit
            end do
            if (umax > 1.02_wp * umaxlo .and. umax > 1.01_wp * ub(jcut)) then
               umax = 0.8_wp * ub(jcut) + 0.2_wp * umax
               if (umax < umaxlo) umax = 1.01_wp * umaxlo
               if (umax < umaxmin) umax = umaxmin
               refine = .true.
            end if
         else if (pstate(nstate1) / pmax > PMIN_UP_QM) then
            umaxlo = umax
            if (1.2_wp * umax < umaxhi) then
               umax = 1.2_wp * umax
               refine = .true.
            else if (umax / umaxhi - 1.0_wp > 0.01_wp) then
               umax = 0.5_wp * (umaxhi + umax)
               refine = .true.
            end if
         end if

         ! Adjust UMIN
         if (pstate(2) / pmax < PMIN_LO_QM) then
            uminlo = umin
            jcut = 1
            do i = 1, nstate1
               jcut = i
               if (pstate(i) / pmax > PMIN_LO_QM) exit
            end do
            if (umin < 0.95_wp * ua(jcut)) then
               umin = 0.2_wp * umin + 0.8_wp * ua(jcut)
               refine = .true.
            end if
         else if (umin > hc .and. pstate(1) / pmax > PMIN_LO_QM) then
            uminhi = umin
            if (0.8_wp * umin > uminlo) then
               umin = max(hc, 0.8_wp * umin)
               refine = .true.
            else if ((umin - uminlo) > 0.01_wp * umin .and. &
                     umin / hc > 20.0_wp) then
               umin = max(0.5_wp * (uminlo + umin), hc)
               refine = .true.
            end if
         end if

         if (.not. refine) then
            converged = .true.
            exit
         end if
         if (iter == MAX_ITER) istat = QM_WINDOW_NOT_CONVERGED
      end do

      if (converged) then
         if (any(pstate /= pstate)) then
            istat = QM_NO_STEADY_STATE
         else
            sum_p = sum(pstate)
            pstate = pstate / sum_p
            call qm_emission(nstate, nlam, isrf_wl_full, cabs_full, &
                             ua, ub, lngu, pstate, t_bins, cur_method, &
                             emission_out)
            ! Energy-conservation gate on the EMISSION SPECTRUM: a correct
            ! steady state emits the power it absorbs.  Test the spectrum
            ! itself (P_abs = int Cabs ISRF dlam, P_em = int lam P_lam dln lam)
            ! and reject a grain that misses by more than 5%.
            block
               real(wp) :: p_abs, p_em, dlnl, wl
               integer  :: ii
               p_abs = 0.0_wp; p_em = 0.0_wp
               do ii = 1, nlam
                  if (ii == 1) then
                     dlnl = 0.5_wp * log(isrf_wl_full(2) / isrf_wl_full(1))
                  else if (ii == nlam) then
                     dlnl = 0.5_wp * log(isrf_wl_full(nlam) / isrf_wl_full(nlam-1))
                  else
                     dlnl = 0.5_wp * log(isrf_wl_full(ii+1) / isrf_wl_full(ii-1))
                  end if
                  wl = isrf_wl_full(ii)
                  p_abs = p_abs + cabs_full(ii) * isrf_full(ii) * wl * dlnl
                  p_em  = p_em  + emission_out(ii) * dlnl
               end do
               if (p_abs > 0.0_wp) ratio_em = p_em / p_abs
               if (any(emission_out /= emission_out)) then
                  istat = QM_NO_STEADY_STATE
               else if (p_abs > 0.0_wp .and. abs(ratio_em - 1.0_wp) > 0.05_wp) then
                  istat = QM_ENERGY_NOT_CONSERVED
               else
                  istat = QM_SOLVED
                  solved = .true.
                  if (present(actual_cooling_out)) actual_cooling_out = cur_method
               end if
            end block
            if (.not. solved) emission_out = 0.0_wp
         end if
      end if
      if (present(status)) status = istat
      if (present(raw_negative_count_out)) raw_negative_count_out = negative_total

      ! Single-grain diagnostic: every outcome, with the reason.
      if (qm_verbose) then
         jpmax = 0
         if (solved) jpmax = maxloc(pstate, dim=1)
         write(0,'(a,es9.2,a,i0,a,i0,a,a,a,l1,a,a,a,f8.5,a,i0,a,i0)') &
            '  QM grain: a[um]=', a_cm * 1.0d4, ' natom=', natom, &
            ' nset=', nset_qm, ' method=', trim(cur_method), ' ok=', solved, &
            ' status=', trim(qm_status_name(istat)), ' Pem/Pabs=', ratio_em, &
            ' iter=', min(iter, MAX_ITER), ' jPk=', jpmax
         write(0,'(a,i0)') '  QM raw negative rates=', negative_total
      end if

      deallocate(isrf_wl_full, isrf_full, cabs_full)
      deallocate(isrf_wl_ds, isrf_ds, cabs_ds)
      deallocate(u, ua, ub, t_bins, lngu)
      deallocate(amatrix, pstate)
   end subroutine qm_solve_grain


   !====================================================================
   ! 19. Name of a qm_solve_grain status code.
   !====================================================================
   function qm_status_name(istat) result(name)
      integer, intent(in) :: istat
      character(len=24) :: name
      select case (istat)
      case (QM_SOLVED);               name = 'solved'
      case (QM_INVALID_INPUT);        name = 'invalid_input'
      case (QM_NO_STEADY_STATE);      name = 'no_steady_state'
      case (QM_WINDOW_NOT_CONVERGED); name = 'window_not_converged'
      case (QM_ENERGY_NOT_CONSERVED); name = 'energy_not_conserved'
      case (QM_ATOM_COUNT_CAP);       name = 'atom_count_cap'
      case (QM_INVALID_RATE);         name = 'invalid_transition_rate'
      case default;                   name = 'unknown'
      end select
   end function qm_status_name

end module stoch_qm_mod
