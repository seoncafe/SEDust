module stoch_qm_mod
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
   !     discrete method, DL01 Section 8.2).
   !   - Intrabin transitions (self-transition correction)
   !   - BiCG sparse linear solver for the transition matrix equation
   !   - Emission from P(state) using the thermal-discrete Planck formula:
   !     sum_J P(J) * B_lam(T_J) * Cabs * 8*pi/lam^4
   !
   ! Two cooling treatments are provided:
   !   'dbdis'  thermal-discrete: all downward transitions f < i retained
   !            (the production reference).
   !   'dbcon'  thermal-continuous: the downward rates are collapsed into
   !            the nearest-neighbour transition, reproducing the
   !            continuous-cooling approximation of the Guhathakurta-Draine
   !            temperature-space recursion.
   !
   ! Units: internally CGS (wavelengths cm, energies erg, cross sections
   ! cm^2, radiation field c*u_lambda erg/cm^3/s). The public API
   ! converts from/to our SI/um convention at the entry/exit points.
   !--------------------------------------------------------------------
   use constants, only: wp
   implicit none
   private

   ! Public API
   public :: qm_solve_grain   ! single-grain QM P(state) solver + emission
   public :: qm_nstate_default
   public :: qm_nisrf_max     ! ISRF downsampling cap for the transition matrix
   public :: qm_method        ! 'dbdis' (thermal-discrete, production) or
                              ! 'dbcon' (thermal-continuous / GD
                              ! nearest-neighbour collapse)
   public :: qm_verbose       ! single-grain stderr diagnostics on/off

   ! Method selector (set once before solving). Default = validated dbdis path.
   character(len=5) :: qm_method = 'dbdis'

   ! Diagnostic-output toggle. Default .true. so the CLI drivers keep the
   ! single-grain stderr lines; dust_emission sets it from the model's
   ! `verbose` field so the library path stays silent by default.
   logical, save :: qm_verbose = .true.

   ! Runtime-settable convention knobs (Draine's production settings: 500 bins, full
   ! 2500-point ISRF grid -- no downsampling).
   integer :: qm_nstate_default = 200  ! enthalpy bins
   integer :: qm_nisrf_max      = 200  ! max ISRF wavelengths for transition matrix

   ! CGS physical constants (matching Draine's published values)
   real(wp), parameter :: HC_CGS  = 6.62607d-27 * 2.99792d10   ! h*c (erg*cm)
   real(wp), parameter :: H_CGS   = 6.62607d-27                ! h (erg*s)
   real(wp), parameter :: C_CGS   = 2.99792d10                 ! c (cm/s)
   real(wp), parameter :: KB_CGS  = 1.38065d-16                ! k_B (erg/K)
   real(wp), parameter :: PI_CGS  = 3.141593d0
   real(wp), parameter :: EV2ERG  = 1.60218d-12                ! eV -> erg

   ! Simpson integration tolerance
   real(wp), parameter :: EPS_SIMP = 1.0d-3

   ! Threshold for sparse matrix storage
   real(wp), parameter :: SPARSE_THRESH = 1.0d-90

contains

   !====================================================================
   ! 1. Parabolic interpolation (3-point Lagrange)
   !    Parabolic interpolation, following Draine's method.
   !====================================================================
   subroutine parab_interp(x, y, n, t, z)
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
   end subroutine parab_interp


   !====================================================================
   ! 2. Interpolate radiation field and cross sections at wavelength
   !    Helper used by all heating/cooling integrands
   !====================================================================
   subroutine interp_rf_cabs(wl_cm, nisrf, isrf_wl, isrf, cabs_arr, &
                             crssct, radfld)
      implicit none
      integer,  intent(in)  :: nisrf
      real(wp), intent(in)  :: wl_cm
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      real(wp), intent(out) :: crssct, radfld

      if (wl_cm >= isrf_wl(1) .and. wl_cm <= isrf_wl(nisrf)) then
         call parab_interp(isrf_wl, cabs_arr, nisrf, wl_cm, crssct)
         call parab_interp(isrf_wl, isrf,     nisrf, wl_cm, radfld)
      else if (wl_cm < isrf_wl(1)) then
         crssct = cabs_arr(1)
         radfld = 0.0_wp
      else
         crssct = cabs_arr(nisrf) * (isrf_wl(nisrf) / wl_cm)**2
         radfld = 0.0_wp
      end if
   end subroutine interp_rf_cabs




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

      wl = HC_CGS / e_ph
      call interp_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      radfld_e = radfld * wl**2 / (H_CGS * C_CGS)

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

      wl = HC_CGS / e_ph
      call interp_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      cradfld_e = radfld * wl**2 / (H_CGS * C_CGS)

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

      wl = HC_CGS / e_ph
      call interp_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      cradfld_e = radfld * wl**2 / (H_CGS * C_CGS)
      hf = (1.0_wp - e_ph / dui) * crssct * cradfld_e
   end function intrabin_heating_funct






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

      elow_cut = HC_CGS / isrf_wl(nisrf)
      eupp_cut = HC_CGS / isrf_wl(1)

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
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) exit
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

      elow_cut = HC_CGS / isrf_wl(nisrf)
      eupp_cut = HC_CGS / isrf_wl(1)

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
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) exit
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

      elow_cut = HC_CGS / isrf_wl(nisrf)
      eupp_cut = HC_CGS / isrf_wl(1)
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
         if (abs(s1) == 0.0_wp .and. abs(s2) == 0.0_wp) exit
         t1 = t2; n = n * 2; hstep = hstep / 2.0_wp; s1 = s2
         if (n > 2**16) exit
      end do
      rate_heating = s2
   end subroutine simpson_intrabin_heating






   !====================================================================
   ! 9b. Thermal cooling integrand (THERM_COOLING_FUNCT)
   !     Follows Draine's method.
   !
   !     The photon-emission integrand divided by the Planck thermal
   !     factor (exp(E/kT) - 1). This is the 'dbdis' (thermal-discrete)
   !     method of DL01 Section 8.2: the downward transition rate is set
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

      boltz_arg = e_ph / (KB_CGS * t_upper)
      if (boltz_arg > 500.0_wp) return   ! exp(-500) ~ 0

      wl = HC_CGS / e_ph
      call interp_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      radfld_e = radfld * wl**2 / (HC_CGS * C_CGS)

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
              (8.0_wp * PI_CGS / (H_CGS * HC_CGS**2)) * &
              (1.0_wp + HC_CGS**3 / (8.0_wp * PI_CGS * e_ph**3) * radfld_e) / &
              (exp(boltz_arg) - 1.0_wp)
      else
         if (e_ph < e1 .or. e_ph >= e4) then
            cf = 0.0_wp
         else
            cf = (e_ph**3) * crssct * &
                 (8.0_wp * PI_CGS / (H_CGS * HC_CGS**2)) * &
                 (1.0_wp + HC_CGS**3 / (8.0_wp * PI_CGS * e_ph**3) * radfld_e) / &
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

      boltz_arg = e_ph / (KB_CGS * t_upper)
      if (boltz_arg > 500.0_wp) return

      wl = HC_CGS / e_ph
      call interp_rf_cabs(wl, nisrf, isrf_wl, isrf, cabs_arr, crssct, radfld)
      radfld_e = radfld * wl**2 / (HC_CGS * C_CGS)

      cf = (8.0_wp * PI_CGS / (H_CGS * HC_CGS**2)) * &
           (1.0_wp - e_ph / dui) * (e_ph**3) * crssct * &
           (1.0_wp + HC_CGS**3 / (8.0_wp * PI_CGS * e_ph**3) * radfld_e) / &
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

      real(wp) :: hstep, t1, t2, s1, s2, p, xval
      integer  :: n, k

      ! Short-circuit if E1 > 20*kT (thermal emission negligible)
      if (e1_io > 20.0_wp * KB_CGS * t_upper) then
         rate_cooling = 0.0_wp; return
      end if

      if (e1_io == e4_io) then
         rate_cooling = 0.0_wp; return
      end if

      n = 1
      hstep = e4_io - e1_io
      t1 = 0.5_wp * hstep * ( &
           therm_cooling_funct(e1_io, t_upper, nisrf, isrf_wl, isrf, cabs_arr, &
                               e1_io, e2_io, e3_io, e4_io, dui, duf) + &
           therm_cooling_funct(e4_io, t_upper, nisrf, isrf_wl, isrf, cabs_arr, &
                               e1_io, e2_io, e3_io, e4_io, dui, duf))
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

      real(wp) :: hstep, t1, t2, s1, s2, p, xval
      integer  :: n, k

      ! Short-circuit if E1 > 20*kT
      if (e1_io > 20.0_wp * KB_CGS * t_upper) then
         rate_cooling = 0.0_wp; return
      end if

      if (e1_io == e2_io) then
         rate_cooling = 0.0_wp; return
      end if

      n = 1
      hstep = e2_io - e1_io
      t1 = 0.5_wp * hstep * ( &
           intrabin_therm_cooling_funct(e1_io, dui, t_upper, nisrf, isrf_wl, isrf, cabs_arr) + &
           intrabin_therm_cooling_funct(e2_io, dui, t_upper, nisrf, isrf_wl, isrf, cabs_arr))
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
                                     delta_u2, nstate, &
                                     nisrf, isrf_wl, isrf, cabs_arr, &
                                     afi_cooling)
      implicit none
      integer,  intent(in)  :: ibin, fbin, nstate, nisrf
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
   !     Transition-matrix construction, 'dbdis' branch, following Draine's method.
   !====================================================================
   subroutine build_transition_matrix(nstate, u, ua, ub, t_bins, &
                                      nisrf, isrf_wl, isrf, cabs_arr, &
                                      method, amatrix, pstate_gd)
      implicit none
      integer,  intent(in)  :: nstate, nisrf
      real(wp), intent(in)  :: u(nstate), ua(nstate), ub(nstate)
      real(wp), intent(in)  :: t_bins(nstate)
      real(wp), intent(in)  :: isrf_wl(nisrf), isrf(nisrf), cabs_arr(nisrf)
      character(len=*), intent(in) :: method   ! 'dbdis' or 'dbcon' (thread-private)
      real(wp), intent(out) :: amatrix(nstate, nstate)
      real(wp), intent(out) :: pstate_gd(nstate) ! Guhathakurta-Draine preconditioner

      integer  :: ibin, fbin
      real(wp) :: afi_val, adiag, elow_lim, eupp_lim
      real(wp) :: delta_u2

      amatrix = 0.0_wp

      ! Pre-compute photon energy limits from the ISRF wavelength range
      elow_lim = HC_CGS / isrf_wl(nisrf)   ! lowest photon energy
      eupp_lim = HC_CGS / isrf_wl(1)       ! highest photon energy

      ! ---- UPWARD (heating) transitions: ibin -> fbin, fbin > ibin ----
      ! These are identical for both methods (dbdis/dbcon).
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
         if (t_bins(ibin) <= 0.0_wp) cycle

         do fbin = ibin - 1, 1, -1
            if (fbin /= 1 .and. ub(fbin) <= ua(fbin)) cycle
            if (u(ibin) <= u(fbin)) cycle

            ! Thermal-discrete (dbdis): Planck-suppressed beyond ~20kT,
            ! so distant downward transitions are negligible and skipped.
            if (ua(ibin) - ub(fbin) > 20.0_wp * KB_CGS * t_bins(ibin)) exit

            call calc_therm_cooling_afi(ibin, t_bins(ibin), &
                                        u(ibin), ua(ibin), ub(ibin), &
                                        fbin, u(fbin), ua(fbin), ub(fbin), &
                                        delta_u2, nstate, &
                                        nisrf, isrf_wl, isrf, cabs_arr, &
                                        afi_val)
            amatrix(fbin, ibin) = afi_val
         end do
      end do

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

      pstate_gd = 0.0_wp
   end subroutine build_transition_matrix


   !====================================================================
   ! 13. Sparse matrix routines + BiCG solver
   !     Biconjugate-gradient sparse solver.
   !====================================================================

   ! Convert dense matrix to sparse row-indexed storage. ok=.false. signals a
   ! storage overflow (more off-diagonal entries than nmax); the caller then
   ! treats the grain as unsolved and falls back to the GD solver.
   subroutine dense_to_sparse(a, n, thresh, nmax, sa, ija, ok)
      implicit none
      integer,  intent(in)  :: n, nmax
      real(wp), intent(in)  :: a(n, n), thresh
      real(wp), intent(out) :: sa(nmax)
      integer,  intent(out) :: ija(nmax)
      logical,  intent(out) :: ok
      integer :: i, j, k

      ok = .true.
      ! Diagonal
      do j = 1, n
         sa(j) = a(j, j)
      end do
      ija(1) = n + 2
      k = n + 1
      do i = 1, n
         do j = 1, n
            if (abs(a(i, j)) >= thresh .and. i /= j) then
               k = k + 1
               if (k > nmax) then
                  if (qm_verbose) &
                     write(*,'(a,i0,a,i0)') 'stoch_qm: sparse overflow k=', k, ' nmax=', nmax
                  ok = .false.
                  return
               end if
               sa(k) = a(i, j)
               ija(k) = j
            end if
         end do
         ija(i + 1) = k + 1
      end do
   end subroutine dense_to_sparse


   ! Sparse matrix * vector: b = A * x
   subroutine sp_matvec(nmax, sa, ija, x, b, n)
      implicit none
      integer,  intent(in)  :: nmax, n
      real(wp), intent(in)  :: sa(nmax), x(n)
      integer,  intent(in)  :: ija(nmax)
      real(wp), intent(out) :: b(n)
      integer :: i, k

      do i = 1, n
         b(i) = sa(i) * x(i)
         do k = ija(i), ija(i+1) - 1
            b(i) = b(i) + sa(k) * x(ija(k))
         end do
      end do
   end subroutine sp_matvec


   ! Sparse transpose matrix * vector: b = A^T * x
   subroutine sp_matvec_t(nmax, sa, ija, x, b, n)
      implicit none
      integer,  intent(in)  :: nmax, n
      real(wp), intent(in)  :: sa(nmax), x(n)
      integer,  intent(in)  :: ija(nmax)
      real(wp), intent(out) :: b(n)
      integer :: i, j, k

      do i = 1, n
         b(i) = sa(i) * x(i)
      end do
      do i = 1, n
         do k = ija(i), ija(i+1) - 1
            j = ija(k)
            b(j) = b(j) + sa(k) * x(i)
         end do
      end do
   end subroutine sp_matvec_t


   ! Diagonal preconditioner solve: x = b / diag(A)
   ! Guarded against zero diagonal elements.
   subroutine sp_diag_solve(nmax, sa, n, b, x)
      implicit none
      integer,  intent(in)  :: nmax, n
      real(wp), intent(in)  :: sa(nmax), b(n)
      real(wp), intent(out) :: x(n)
      integer :: i
      do i = 1, n
         if (abs(sa(i)) > 0.0_wp) then
            x(i) = b(i) / sa(i)
         else
            x(i) = 0.0_wp
         end if
      end do
   end subroutine sp_diag_solve


   ! L2 norm
   function vec_norm(n, sx) result(snrm)
      implicit none
      integer,  intent(in) :: n
      real(wp), intent(in) :: sx(n)
      real(wp) :: snrm
      snrm = sqrt(sum(sx**2))
   end function vec_norm


   ! BiCG solver for sparse system A*x = b
   ! LINBCG: the Numerical Recipes biconjugate-gradient solver
   subroutine linbcg(nmax, sa, ija, n, b, x, tol_in, itmax, iter_out, err_out)
      implicit none
      integer,  intent(in)    :: nmax, n, itmax
      real(wp), intent(in)    :: sa(nmax), b(n), tol_in
      integer,  intent(in)    :: ija(nmax)
      real(wp), intent(inout) :: x(n)
      integer,  intent(out)   :: iter_out
      real(wp), intent(out)   :: err_out

      real(wp), allocatable :: p(:), pp(:), r(:), rr(:), z(:), zz(:)
      real(wp) :: ak, akden, bk, bkden, bknum, bnrm, &
                  zm1nrm, znrm
      real(wp), parameter :: EPS_BCG = 1.0d-14
      integer :: j, iter

      allocate(p(n), pp(n), r(n), rr(n), z(n), zz(n))

      iter_out = 0

      ! r = b - A*x
      call sp_matvec(nmax, sa, ija, x, r, n)
      do j = 1, n
         r(j) = b(j) - r(j)
         rr(j) = r(j)
      end do

      znrm = 1.0_wp
      ! ITOL=1: use ||b|| as normalization
      bnrm = vec_norm(n, b)
      if (bnrm == 0.0_wp) bnrm = 1.0_wp

      call sp_diag_solve(nmax, sa, n, r, z)

      bkden = 1.0_wp
      do iter = 1, itmax
         iter_out = iter
         zm1nrm = znrm
         call sp_diag_solve(nmax, sa, n, rr, zz)
         bknum = 0.0_wp
         do j = 1, n
            bknum = bknum + z(j) * rr(j)
         end do
         if (iter == 1) then
            do j = 1, n
               p(j)  = z(j)
               pp(j) = zz(j)
            end do
         else
            if (abs(bkden) > 0.0_wp) then
               bk = bknum / bkden
            else
               bk = 0.0_wp
            end if
            do j = 1, n
               p(j)  = bk * p(j)  + z(j)
               pp(j) = bk * pp(j) + zz(j)
            end do
         end if
         bkden = bknum
         call sp_matvec(nmax, sa, ija, p, z, n)
         akden = 0.0_wp
         do j = 1, n
            akden = akden + z(j) * pp(j)
         end do
         if (abs(akden) > 0.0_wp) then
            ak = bknum / akden
         else
            ak = 0.0_wp
         end if
         call sp_matvec_t(nmax, sa, ija, pp, zz, n)
         do j = 1, n
            x(j)  = x(j)  + ak * p(j)
            r(j)  = r(j)  - ak * z(j)
            rr(j) = rr(j) - ak * zz(j)
         end do
         call sp_diag_solve(nmax, sa, n, r, z)

         err_out = vec_norm(n, r) / bnrm
         if (err_out <= tol_in) exit
      end do

      deallocate(p, pp, r, rr, z, zz)
   end subroutine linbcg




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

      ! Convert Debye temperatures to cm^-1 (T / (hc/k) = T / 1.4388 cm^-1)
      em2d = TD2D / 1.4388_wp
      em3d = TD3D / 1.4388_wp

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
   ! 14g. T_FROM_U_BISECT -- invert U(T) -> T for enthalpy bin labeling.
   !      Uses our Debye enthalpy_DL01 via bisection.
   !      U_target is in erg; returned T is in K.
   !====================================================================
   function T_from_U_bisect(U_target, grain_type, natom, a_cm) result(T)
      use enthalpy, only: enthalpy_DL01
      implicit none
      real(wp),         intent(in) :: U_target, a_cm
      character(len=*), intent(in) :: grain_type
      integer,          intent(in) :: natom
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
   !        nset_out    = number of mode-tracking bins used
   !====================================================================
   subroutine build_enthalpy_bins_qm(grain_type, a_cm, natom, nc, nh, &
                                      umin, umax, nstate, &
                                      u, ua, ub, t_bins, nset_out)
      implicit none
      character(len=*), intent(in)  :: grain_type
      real(wp),         intent(in)  :: a_cm, umin, umax
      integer,          intent(in)  :: natom, nc, nh, nstate
      real(wp),         intent(out) :: u(nstate), ua(nstate), ub(nstate)
      real(wp),         intent(out) :: t_bins(nstate)
      integer,          intent(out) :: nset_out

      integer  :: nmodes, i, nset, nset_max, jset
      real(wp) :: de, dlgu, umax_cm1
      real(wp), allocatable :: emodes(:), emodes_new(:)
      logical  :: skip_modes

      ! ---- Step 1: Get vibrational mode spectrum ----
      ! With UMIN /= 0 the bins are purely log-spaced, so the explicit mode
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

      ! ---- Step 2: Determine DE (bin-boundary snapping quantum) ----
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
            ! Fallback: too few modes for mode-tracking; use log-spaced
            nset = 0; nset_out = 0
            call fallback_log_bins(nstate, umin, umax, u, ua, ub)
            deallocate(emodes_new)
            goto 500   ! jump to the T(U) computation
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
         umax_cm1 = umax / HC_CGS
         jset = nstate   ! fallback

         do i = nset + 1, nstate
            u(i) = u(nset) + (ub(nset) - ua(nset)) * real(i - nset, wp)
            dlgu = log10(umax_cm1 / u(i)) / real(nstate - i, wp)
            if (nstate - i > 0 .and. &
                (u(i) * (10.0_wp**dlgu) - u(i)) > (u(i) - u(i - 1))) then
               jset = i
               exit
            end if
            if (i == nstate) then
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
            u(i)  = u(i)  * HC_CGS
            ua(i) = ua(i) * HC_CGS
            ub(i) = ub(i) * HC_CGS
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
      ! Use the Debye model, for consistency with the GD solver's
      ! enthalpy_DL01.
500   continue
      do i = 1, nstate
         t_bins(i) = T_from_U_bisect(u(i), grain_type, natom, a_cm)
      end do

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


      ! Fallback: simple log-spaced bins when mode-tracking is infeasible
      subroutine fallback_log_bins(ns, umin_f, umax_f, u_f, ua_f, ub_f)
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
      end subroutine fallback_log_bins



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

      select case (trim(grain_type))
      case ('sil')
         natom = nint((4.0_wp * PI_CGS / 3.0_wp) * a_cm**3 * &
                      3.5_wp * 7.0_wp / (172.0_wp * 1.66d-24))
      case ('pah')
         nc_real = 472.0d21 * a_cm**3
         natom = max(nint(nc_real), 3)
      case default
         natom = nint((4.0_wp * PI_CGS / 3.0_wp) * a_cm**3 * &
                      3.5_wp * 7.0_wp / (172.0_wp * 1.66d-24))
      end select
      natom = max(natom, 3)
   end function compute_natom


   !====================================================================
   ! 17. QM emission from P(state) using the thermal-discrete formula.
   !
   !       EMISSION(I) = sum_J P(J) * (8*pi*hcc/lambda^4) * Cabs
   !                     / (exp(hc/(kT_J * lambda)) - 1)
   !     for bins J where h*nu <= UB(J) (photon energy < bin energy).
   !
   !     This is the thermal-discrete approximation of DL01 Section 6.1:
   !     each bin emits as a blackbody at its temperature T(J), which is
   !     the temperature corresponding to the bin's representative energy
   !     via the T(E) relation of DL01 eq. (32).
   !
   !     NOTE: The (1+J/B0) stimulated-emission factor is omitted
   !     (HD23 publishes net emission only).
   !====================================================================
   subroutine qm_emission(nstate, nisrf, isrf_wl_cm, cabs_cm2, &
                           u, ua, ub, pstate, t_bins, &
                           method, emission)
      implicit none
      integer,  intent(in)  :: nstate, nisrf
      real(wp), intent(in)  :: isrf_wl_cm(nisrf), cabs_cm2(nisrf)
      real(wp), intent(in)  :: u(nstate), ua(nstate), ub(nstate)
      real(wp), intent(in)  :: pstate(nstate)
      real(wp), intent(in)  :: t_bins(nstate)
      character(len=*), intent(in) :: method   ! 'dbdis' or 'dbcon' (thread-private)
      real(wp), intent(out) :: emission(nisrf)

      integer  :: i, j
      real(wp) :: hnu, term, hcc_val, dj

      hcc_val = H_CGS * C_CGS * C_CGS   ! h * c^2

      emission = 0.0_wp

      do i = 1, nisrf
         hnu = HC_CGS / isrf_wl_cm(i)

         ! Thermal-discrete (dbdis) kernel: blackbody at bin temperature T_j.
         do j = 2, nstate
            if (pstate(j) <= 0.0_wp) cycle
            if (t_bins(j) <= 0.0_wp) cycle
            ! Photon energy must not exceed the upper bound of this bin
            if (hnu > ub(j)) cycle

            term = HC_CGS / (isrf_wl_cm(i) * KB_CGS * t_bins(j))
            if (term < 500.0_wp) then
               dj = (8.0_wp * PI_CGS * hcc_val / isrf_wl_cm(i)**4) * &
                    cabs_cm2(i) * pstate(j) / (exp(term) - 1.0_wp)
               emission(i) = emission(i) + dj
            end if
         end do
      end do
   end subroutine qm_emission


   !====================================================================
   ! 18. Top-level single-grain QM driver
   !
   ! Given a grain with:
   !   - cross section Cabs(lambda) in cm^2
   !   - radiation field J_lam in SI W/m^3/sr at wavelengths lam_um
   !   - enthalpy table H_wide(T) in erg and temperature grid T_wide(K)
   !   - equilibrium enthalpy EEQ (erg)
   !   - grain type ('sil' or 'pah')
   !   - grain radius a_cm
   !
   ! Returns emission(NLAM) in CGS erg/s/cm per grain, matching the
   ! unit convention of Draine's emission array.
   !
   ! The caller (sed_grain_loop) accumulates:
   !   Jout += dn_pop(ir) * emission_qm * (units factor)
   !====================================================================
   subroutine qm_solve_grain(nlam, lam_um, cabs_cm2, j_lam_si, &
                              nt_wide, t_wide, h_wide, &
                              teq, eeq, eeqss, &
                              a_cm, grain_type, &
                              emission_out, solved)
      use sed_mathlib, only: interp
      implicit none
      integer,          intent(in)  :: nlam, nt_wide
      real(wp),         intent(in)  :: lam_um(nlam), cabs_cm2(nlam)
      real(wp),         intent(in)  :: j_lam_si(nlam)
      real(wp),         intent(in)  :: t_wide(nt_wide), h_wide(nt_wide)
      real(wp),         intent(in)  :: teq, eeq, eeqss
      real(wp),         intent(in)  :: a_cm
      character(len=*), intent(in)  :: grain_type
      real(wp),         intent(out) :: emission_out(nlam)
      logical,          intent(out) :: solved

      ! Local variables
      integer :: nstate, nstate1, nmax, natom
      integer :: i, j, fbin, iter, iter_total, jcut, jpmax
      integer :: nc_pah, nh_pah, nset_qm
      real(wp) :: umin, umax, umaxmin, umaxhi, umaxlo, uminhi, uminlo
      real(wp) :: pmax, term, sum_p, err_bcg, tol_bcg
      logical  :: refine, bicg_ok, sparse_ok
      integer  :: itmax_bcg, n_retry
      integer, parameter :: MAX_BCG_RETRY = 30
      ! Full-resolution CGS arrays for emission computation
      real(wp), allocatable :: isrf_wl_full(:), isrf_full(:), cabs_full(:)
      ! Downsampled CGS arrays for transition matrix (much faster)
      integer :: nisrf_ds
      real(wp), allocatable :: isrf_wl_ds(:), isrf_ds(:), cabs_ds(:)

      real(wp), allocatable :: u(:), ua(:), ub(:), t_bins(:)
      real(wp), allocatable :: amatrix(:,:), pstate_gd(:), pstate(:)
      real(wp), allocatable :: amatrix1(:,:), pstate1(:), rhs(:)
      real(wp), allocatable :: sa(:)
      integer,  allocatable :: ija(:)

      real(wp), parameter :: PMIN_LO_QM = 1.0d-13
      real(wp), parameter :: PMIN_UP_QM = 1.0d-13
      integer,  parameter :: MAX_ITER = 10

      ! Largest atom count for which the explicit-mode treatment (PAH_MODES /
      ! sil modes, ~3*natom modes) is feasible. Above this we return
      ! not-solved so the caller's GD fallback (Debye H(T), no mode array)
      ! handles the grain -- such large grains sit near the equilibrium
      ! boundary anyway. Prevents the ~3e7-element allocation seen for the
      ! DL07 carbonaceous population, which extends to large a with non-zero dn.
      integer, parameter :: NATOM_QM_MAX = 1000000

      solved = .false.
      emission_out = 0.0_wp

      ! Determine number of enthalpy bins
      nstate = qm_nstate_default

      ! Number of atoms in the grain
      natom = compute_natom(a_cm, grain_type)

      ! Compute NC, NH for PAH grains (needed by mode spectrum builder)
      nc_pah = 0; nh_pah = 0
      if (trim(grain_type) == 'pah') then
         call pah_size_atoms(a_cm, nc_pah, nh_pah)
      end if

      ! Feasibility cap: explicit-mode construction (~3*natom modes) is only
      ! needed when UMIN = 0 (mode-aligned first bins; eeq < 0.1*eeqss).  With
      ! UMIN = eeq/5 > 0 the bins are log-spaced and modes are skipped, so
      ! large grains are feasible.
      if (max(natom, nc_pah) > NATOM_QM_MAX .and. eeq < 0.1_wp * eeqss) then
         solved = .false.
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
         isrf_full(i)    = 4.0_wp * PI_CGS * 10.0_wp * j_lam_si(i)
         cabs_full(i)    = cabs_cm2(i)
      end do

      ! Downsample to NISRF_MAX_QM points (log-spaced) for transition matrix
      ! This is the key optimization: Draine's method uses ~200 wavelengths, not 1129.
      if (nlam > qm_nisrf_max) then
         nisrf_ds = qm_nisrf_max
         allocate(isrf_wl_ds(nisrf_ds), isrf_ds(nisrf_ds), cabs_ds(nisrf_ds))
         do i = 1, nisrf_ds
            ! Map to index in full array via equal log-spacing
            j = 1 + nint(real(i - 1, wp) / real(nisrf_ds - 1, wp) * real(nlam - 1, wp))
            j = max(1, min(j, nlam))
            isrf_wl_ds(i) = isrf_wl_full(j)
            isrf_ds(i)    = isrf_full(j)
            cabs_ds(i)    = cabs_full(j)
         end do
      else
         nisrf_ds = nlam
         allocate(isrf_wl_ds(nisrf_ds), isrf_ds(nisrf_ds), cabs_ds(nisrf_ds))
         isrf_wl_ds = isrf_wl_full
         isrf_ds    = isrf_full
         cabs_ds    = cabs_full
      end if

      ! Initial UMIN/UMAX guesses (following Draine's method)
      umaxmin = 13.65_wp * EV2ERG
      umax = max(13.6_wp * EV2ERG + 2.0_wp * eeq, umaxmin)
      if (eeq < 0.1_wp * eeqss) then
         umin = 0.0_wp
      else
         umin = eeq / 5.0_wp
      end if
      umaxhi = 1.0d70 * HC_CGS
      umaxlo = 0.0_wp
      uminhi = 1.0d70 * HC_CGS
      uminlo = 0.0_wp

      ! Allocate state arrays
      allocate(u(nstate), ua(nstate), ub(nstate), t_bins(nstate))
      allocate(amatrix(nstate, nstate), pstate_gd(nstate), pstate(nstate))

      nstate1 = nstate - 1
      nmax = (nstate1 + 1)**2 + 1
      allocate(amatrix1(nstate1, nstate1), pstate1(nstate1), rhs(nstate1))
      allocate(sa(nmax), ija(nmax))

      ! Iterative UMIN/UMAX refinement loop (following Draine's method)
      do iter = 1, MAX_ITER

         ! Build enthalpy bins with mode-aware QM construction
         call build_enthalpy_bins_qm(grain_type, a_cm, natom, nc_pah, nh_pah, &
                                     umin, umax, nstate, &
                                     u, ua, ub, t_bins, nset_qm)

         ! Build transition matrix (using downsampled ISRF for speed)
         call build_transition_matrix(nstate, u, ua, ub, t_bins, &
                                      nisrf_ds, isrf_wl_ds, isrf_ds, cabs_ds, &
                                      qm_method, amatrix, pstate_gd)

         ! Check for low-U limit: direct solve without BiCG
         fbin = 1
         do i = 1, nstate
            if (u(i) < 2.0d-12) fbin = i   ! 2e-12 erg ~ 12.5 eV
         end do
         if (amatrix(fbin,fbin) /= 0.0_wp) then
            term = amatrix(1,1) / amatrix(fbin,fbin)
         else
            term = 1.0_wp
         end if

         if (abs(term) < 1.0d-8) then
            ! Low-U direct solve
            pstate(1) = 1.0_wp
            if (amatrix(nstate,nstate) /= 0.0_wp) then
               pstate(nstate) = -amatrix(nstate,1) / amatrix(nstate,nstate)
            else
               pstate(nstate) = 0.0_wp
            end if
            do fbin = nstate - 1, 2, -1
               term = amatrix(fbin, 1)
               do i = fbin + 1, nstate
                  term = term + pstate(i) * amatrix(fbin, i)
               end do
               if (amatrix(fbin,fbin) /= 0.0_wp) then
                  pstate(fbin) = -term / amatrix(fbin, fbin)
               else
                  pstate(fbin) = 0.0_wp
               end if
            end do
            ! Normalize
            term = sum(pstate)
            if (term > 0.0_wp) pstate = pstate / term
            pstate1(1:nstate1) = pstate(1:nstate1)
            solved = .true.
            exit   ! No UMIN/UMAX refinement needed for low-U

         else
            ! Full BiCG solve
            ! Reduce system to (NSTATE-1) equations: P_N = 1 - sum(P_1..N-1)
            do fbin = 1, nstate1
               rhs(fbin) = -amatrix(fbin, nstate)
               pstate1(fbin) = 0.0_wp
               do i = 1, nstate1
                  amatrix1(fbin, i) = amatrix(fbin, i) - amatrix(fbin, nstate)
               end do
            end do

            ! Sparse storage
            sa = 0.0_wp
            ija = 0
            call dense_to_sparse(amatrix1, nstate1, SPARSE_THRESH, nmax, sa, ija, sparse_ok)
            if (.not. sparse_ok) then
               ! storage overflow: zero the state so this grain is reported
               ! unsolved (caller falls back to GD), same as a BiCG failure.
               pstate1(1:nstate1) = 0.0_wp
               pstate = 0.0_wp
               exit
            end if

            ! BiCG solve with convergence enforcement:
            ! grow ITMAX until
            ! err <= tol (warm start), then verify the row residuals of
            ! the linear system; rerun while any row with P > 1e-10 has
            ! |residual|/|A(1,1)| > 1e-6. The original port accepted the
            ! first linbcg result unchecked, which let stagnated solves
            ! (garbage P, factors of 1e4-1e7) through at NSTATE=500 or
            ! with the full 2500-point ISRF grid.
            tol_bcg = 1.0d-15
            itmax_bcg = 1000
            bicg_ok = .false.
            do n_retry = 1, MAX_BCG_RETRY
               itmax_bcg = itmax_bcg + 100
               call linbcg(nmax, sa, ija, nstate1, rhs, pstate1, &
                           tol_bcg, itmax_bcg, iter_total, err_bcg)
               if (err_bcg > tol_bcg) cycle
               bicg_ok = .true.
               do fbin = 1, nstate1
                  if (pstate1(fbin) > 1.0d-10) then
                     term = -rhs(fbin)
                     do i = 1, nstate1
                        term = term + amatrix1(fbin, i) * pstate1(i)
                     end do
                     term = abs(term / amatrix1(1, 1))
                     if (term > 1.0d-6) then
                        bicg_ok = .false.
                        exit
                     end if
                  end if
               end do
               if (bicg_ok) exit
            end do
            if (.not. bicg_ok) then
               ! persistent BiCG failure: zero the state so this grain is
               ! reported unsolved (caller falls back to GD).
               pstate1(1:nstate1) = 0.0_wp
               pstate = 0.0_wp
               exit
            end if

            ! Clamp negative values
            do i = 1, nstate1
               if (pstate1(i) < 0.0_wp) pstate1(i) = 0.0_wp
            end do

            sum_p = sum(pstate1(1:nstate1))

            ! Reconstruct full PSTATE including bin NSTATE
            pstate(1:nstate1) = pstate1(1:nstate1)
            pstate(nstate) = 0.0_wp
            do i = 1, nstate1
               pstate(nstate) = pstate(nstate) + amatrix(nstate, i) * pstate1(i)
            end do
            if (amatrix(nstate,nstate) /= 0.0_wp) then
               pstate(nstate) = -pstate(nstate) / amatrix(nstate, nstate)
            else
               pstate(nstate) = 0.0_wp
            end if
            if (pstate(nstate) < 0.0_wp) pstate(nstate) = 0.0_wp

            ! Find peak
            pmax = maxval(pstate1(1:nstate1))
            if (pmax <= 0.0_wp) exit   ! degenerate

            jpmax = maxloc(pstate1(1:nstate1), dim=1)

            ! Adjust UMAX
            refine = .false.

            if (pstate1(nstate1) / pmax <= PMIN_UP_QM .and. umax > umaxmin) then
               umaxhi = umax
               jcut = nstate1
               do i = nstate1, 1, -1
                  jcut = i
                  if (ub(i) < umaxmin) exit
                  if (pstate1(i) / pmax > PMIN_UP_QM) exit
               end do
               if (umax > 1.02_wp * umaxlo .and. umax > 1.01_wp * ub(jcut)) then
                  umax = 0.8_wp * ub(jcut) + 0.2_wp * umax
                  if (umax < umaxlo) umax = 1.01_wp * umaxlo
                  if (umax < umaxmin) umax = umaxmin
                  refine = .true.
               end if
            else if (pstate1(nstate1) / pmax > PMIN_UP_QM) then
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
            if (nstate1 >= 2 .and. pstate1(2) / pmax < PMIN_LO_QM) then
               uminlo = umin
               jcut = 1
               do i = 1, nstate1
                  jcut = i
                  if (pstate1(i) / pmax > PMIN_LO_QM) exit
               end do
               if (umin < 0.95_wp * ua(jcut)) then
                  umin = 0.2_wp * umin + 0.8_wp * ua(jcut)
                  refine = .true.
               end if
            else if (umin > HC_CGS .and. pstate1(1) / pmax > PMIN_LO_QM) then
               uminhi = umin
               if (0.8_wp * umin > uminlo) then
                  umin = max(HC_CGS, 0.8_wp * umin)
                  refine = .true.
               else if ((umin - uminlo) > 0.01_wp * umin .and. &
                        umin / HC_CGS > 20.0_wp) then
                  umin = max(0.5_wp * (uminlo + umin), HC_CGS)
                  refine = .true.
               end if
            end if

            if (.not. refine) then
               solved = .true.
               exit
            end if
         end if
      end do

      ! If iteration did not converge, accept last result if pmax > 0
      if (.not. solved) then
         pmax = maxval(abs(pstate))
         if (pmax > 0.0_wp) solved = .true.
      end if

      ! Sanity check: if pstate contains NaN, fall back to the
      ! Guhathakurta-Draine preconditioner estimate
      if (solved) then
         do i = 1, nstate
            if (pstate(i) /= pstate(i)) then
               ! Use the GD estimate as fallback
               pstate = pstate_gd
               exit
            end if
         end do
         ! Double-check the fallback for NaN too
         do i = 1, nstate
            if (pstate(i) /= pstate(i)) then
               solved = .false.
               exit
            end if
         end do
      end if

      ! Renormalize pstate so sum(pstate) = 1 (may have drifted from
      ! clamping negatives and from finite BiCG tolerance).
      if (solved) then
         sum_p = sum(pstate(1:nstate))
         if (sum_p > 0.0_wp) pstate = pstate / sum_p
      end if

      ! Compute emission spectrum using QM formula (full-resolution ISRF)
      if (solved) then
         call qm_emission(nstate, nlam, isrf_wl_full, cabs_full, &
                          u, ua, ub, pstate, t_bins, qm_method, emission_out)
         ! Replace any NaN in emission with zero
         do i = 1, nlam
            if (emission_out(i) /= emission_out(i)) emission_out(i) = 0.0_wp
         end do
      end if

      ! Energy-conservation gate on the EMISSION SPECTRUM (ISRFabs vs
      ! P_em). A correct
      ! steady state emits exactly the power it absorbs; a stagnated BiCG
      ! solution that passes the row-residual test can still produce a
      ! spurious emission spectrum (e.g. one silicate grain at nstate=300
      ! emitting ~40x in the MIR). The matrix-transition energy balances
      ! for such a solution -- so we test the actual spectrum directly:
      !   P_abs = integral Cabs(lam) * ISRF(lam) dlam
      !   P_em  = integral (nu*P_nu) dln(lam)
      ! and fall back to GD if they disagree by > 5%.
      if (solved) then
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
            if (p_abs > 0.0_wp) then
               if (abs(p_em / p_abs - 1.0_wp) > 0.05_wp) then
                  solved = .false.
                  emission_out = 0.0_wp
               end if
            end if
         end block
      end if

      ! Single-grain diagnostic
      if (qm_verbose) then
         if (solved) then
            ! Dump transition matrix diagnostics for first small grain
            block
               real(wp) :: sum_heat_1, cool_21, p1, p2, pmax_diag
               integer  :: jpmax_diag

               ! Total heating rate out of bin 1
               sum_heat_1 = 0.0_wp
               do i = 2, nstate
                  sum_heat_1 = sum_heat_1 + amatrix(i, 1)
               end do
               ! Cooling rate from bin 2 to bin 1
               cool_21 = amatrix(1, 2)
               p1 = pstate(1)
               p2 = pstate(2)
               pmax_diag = maxval(pstate(1:nstate))
               jpmax_diag = maxloc(pstate(1:nstate), dim=1)

               write(0,'(a,es9.2,a,i5,a,i3,a,l1,a,i3)') &
                  '  QM grain: a[um]=', a_cm * 1.0d4, &
                  ' natom=', natom, ' nset=', nset_qm, &
                  ' ok=', solved, ' jPk=', jpmax_diag
            end block
         else
            write(0,'(a,es9.2,a,i5,a,i3,a,l1)') &
               '  QM grain: a[um]=', a_cm * 1.0d4, &
               ' natom=', natom, ' nset=', nset_qm, &
               ' ok=', solved
         end if
      end if

      deallocate(isrf_wl_full, isrf_full, cabs_full)
      deallocate(isrf_wl_ds, isrf_ds, cabs_ds)
      deallocate(u, ua, ub, t_bins)
      deallocate(amatrix, pstate_gd, pstate)
      deallocate(amatrix1, pstate1, rhs, sa, ija)
   end subroutine qm_solve_grain

end module stoch_qm_mod
