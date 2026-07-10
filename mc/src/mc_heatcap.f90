module mc_heatcap
   ! Selectable heat-capacity / enthalpy back-ends for the MC P(T) simulation.
   !
   !   'gra_dl01' (default) -- Draine & Li 2001 graphite (Car0 in enthalpy_DL01)
   !                           Two-mode Debye_2(863K) + 2 Debye_2(2504K), plus
   !                           CH stretch modes for small grains (natom < 1e4).
   !
   !   'gra_d85'            -- Draine 1985 original: single Debye(Theta=420 K)
   !                           on N atoms with n = 1.14e23 cm^-3.  Eq. (2.5) of
   !                           Draine & Anderson 1985, ApJ 292, 494.
   !                           Pure Debye-3 internal energy; no CH modes.
   !
   !   'sil_dl01'           -- DL01 silicate ('Sil' in enthalpy_DL01)
   !
   ! Two routines per composition:
   !   U_of_T(T, a_um, comp)   --> internal energy [erg]   (>= 0, U(0)=0)
   !   C_of_T(T, a_um, comp)   --> heat capacity   [erg/K] (>  0)
   !
   ! U is integrated from T=0; C = dU/dT.  For DL01 we wrap enthalpy_DL01
   ! (already integrated from 0); for D85 we evaluate the Debye integral
   ! analytically via the standard Debye_3 function.

   use constants,             only: wp, pi
   use enthalpy,              only: enthalpy_DL01
   use enthalpy_astrodust_mod, only: enthalpy_S1, enthalpy_S2
   implicit none
   private
   public :: U_of_T, C_of_T
   public :: HC_DEFAULT
   character(len=*), parameter :: HC_DEFAULT = 'gra_dl01'

   real(wp), parameter :: kB_cgs = 1.3806488e-16_wp   ! erg / K
   ! Draine 1985 graphite: Debye temperature and atom number density
   real(wp), parameter :: THETA_D85 = 420.0_wp        ! K
   real(wp), parameter :: NDENS_D85 = 1.14e23_wp      ! atoms / cm^3

contains

   function U_of_T(T, a_um, comp) result(U)
      ! Internal vibrational energy [erg] for one grain at temperature T.
      real(wp),         intent(in) :: T, a_um
      character(len=*), intent(in) :: comp
      real(wp) :: U
      real(wp) :: natom, theta_over_T, D3
      if (T <= 0.0_wp) then
         U = 0.0_wp
         return
      end if
      select case (trim(comp))
      case ('gra_dl01', 'pah')
         ! Graphite (DL01 Car0).  Also used for PAH per HD23 §3.2.
         U = enthalpy_DL01(T, a_um, 'Car0')
      case ('sil_dl01')
         U = enthalpy_DL01(T, a_um, 'Sil ')
      case ('ad_s1_c1', 'ad_s1_c2')
         ! enthalpy_S1 is now silicate-only (no C1/C2 charge split); both
         ! legacy keys map to it. Signature is (T, radius_um).
         U = enthalpy_S1(T, a_um)
      case ('ad_s2')
         U = enthalpy_S2(T, a_um)
      case ('gra_d85')
         natom        = NDENS_D85 * (4.0_wp/3.0_wp) * pi * (a_um * 1.0e-4_wp)**3
         theta_over_T = THETA_D85 / T
         D3           = debye3_function(theta_over_T)
         ! Internal energy of a Debye solid: U = 3 N k_B T * D3(Theta/T)
         U = 3.0_wp * natom * kB_cgs * T * D3
      case default
         write(*,'(a,a)') 'mc_heatcap::U_of_T: unknown comp ', trim(comp)
         stop 1
      end select
   end function U_of_T


   function C_of_T(T, a_um, comp) result(C)
      ! Heat capacity dU/dT [erg/K].  Uses a small finite difference on
      ! U_of_T -- robust across all back-ends.  T in K, a in um.
      real(wp),         intent(in) :: T, a_um
      character(len=*), intent(in) :: comp
      real(wp) :: C, dT, Tp, Tm
      if (T <= 0.0_wp) then
         C = 0.0_wp
         return
      end if
      dT = max(1.0e-3_wp, T * 1.0e-3_wp)
      Tp = T + dT
      Tm = max(T - dT, 1.0e-6_wp)
      C  = (U_of_T(Tp, a_um, comp) - U_of_T(Tm, a_um, comp)) / (Tp - Tm)
      if (C <= 0.0_wp) C = tiny(1.0_wp)
   end function C_of_T


   pure function debye3_function(x) result(D3)
      ! Debye function of order 3:  D3(x) = (3/x^3) * Int_0^x t^3/(e^t-1) dt
      ! Numerical integration via Simpson's rule on an x-dependent grid.
      ! Accurate to ~1e-6 for x in [1e-3, 50]; asymptotes handled explicitly.
      real(wp), intent(in) :: x
      real(wp) :: D3
      integer  :: i, n
      real(wp) :: h, s, t, f0, fn, integ
      real(wp), parameter :: pi4_15 = 6.493939402266829149_wp   ! pi^4/15
      if (x <= 0.0_wp) then
         D3 = 1.0_wp
         return
      end if
      if (x < 1.0e-3_wp) then
         ! Series: D3(x) -> 1 - 3x/8 + x^2/20 ...
         D3 = 1.0_wp - 0.375_wp*x + x*x/20.0_wp
         return
      end if
      if (x > 25.0_wp) then
         ! Large x: integrand decays exponentially; integral -> pi^4/15
         D3 = 3.0_wp * pi4_15 / x**3
         return
      end if
      n = 400
      if (n/2*2 /= n) n = n + 1   ! ensure even
      h = x / real(n, wp)
      f0 = 0.0_wp                  ! integrand at t=0: lim = 0
      fn = x**3 / (exp(x) - 1.0_wp)
      s = f0 + fn
      do i = 1, n-1
         t = i * h
         if (i == 2*(i/2)) then    ! even index -> factor 2
            s = s + 2.0_wp * t**3 / (exp(t) - 1.0_wp)
         else                       ! odd  index -> factor 4
            s = s + 4.0_wp * t**3 / (exp(t) - 1.0_wp)
         end if
      end do
      integ = s * h / 3.0_wp
      D3 = 3.0_wp / x**3 * integ
   end function debye3_function

end module mc_heatcap
