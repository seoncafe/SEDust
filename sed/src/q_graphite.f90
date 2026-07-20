module q_graphite_mod
   ! Random-orientation-averaged graphite Q_abs(a, lambda).
   !
   ! Implements DL07 eq. 5-7 / HD23 eq. 15-16 graphite branch: the
   ! C_abs(graphite, a) that the PAH transition function ξ_gra(a) blends
   ! into the PAH cross section for a > 50 Å.
   !
   ! Reads the Draine (2003) graphite dielectric functions and computes
   ! the Car0/Car1 cross sections. Uses the D03 graphite
   ! refractive-index tables index_CpaD03 (E ∥ c) and index_CpeD03
   ! (E ⊥ c), with the size-dependent free-electron contribution
   ! computed analytically per radius. Random-orientation average is
   ! Q_abs = (1/3) Q_∥ + (2/3) Q_⊥.

   use constants, only: wp
   use sed_mathlib,   only: interp
   use mie_mod,   only: mie
   implicit none
   private
   public :: q_graphite_abs
   public :: q_graphite_full

   character(len=*), parameter :: F_CPA = '../data/dielectric/index_CpaD03'
   character(len=*), parameter :: F_CPE = '../data/dielectric/index_CpeD03'

   integer,  parameter :: NCPA = 387
   integer,  parameter :: NCPE = 384
   real(wp), parameter :: T_GRAIN = 20.0_wp           ! K, Draine default
   real(wp), parameter :: PI_LOC  = 3.141592653589793238462643383279502884197_wp

   ! Base dielectric tables (without free-electron contribution)
   logical  :: loaded = .false.
   real(wp) :: cpa_eV(NCPA), cpa_eps1(NCPA), cpa_eps2(NCPA), cpa_wavl(NCPA)
   real(wp) :: cpe_eV(NCPE), cpe_eps1(NCPE), cpe_eps2(NCPE), cpe_wavl(NCPE)

   ! Cached (n, k) tables for current radius. The free-electron
   ! contribution depends on grain size, so we rebuild these only when
   ! `agrain` changes. sed_init walks the radius grid once, so the cache
   ! hit rate is ~99% (one full rebuild per radius, then all wavelengths
   ! reuse the cache).
   real(wp) :: cached_a = -1.0_wp
   real(wp) :: cpa_n(NCPA), cpa_k(NCPA)
   real(wp) :: cpe_n(NCPE), cpe_k(NCPE)

contains

   subroutine load_tables()
      integer  :: i, u
      real(wp) :: ener, rn1, rk, e1, e2

      open(newunit=u, file=F_CPA, status='old', action='read')
      read(u, '(/)')
      do i = 1, NCPA
         read(u, *) ener, rn1, rk, e1, e2
         cpa_eV(i)   = ener
         cpa_eps1(i) = e1
         cpa_eps2(i) = e2
         cpa_wavl(i) = 1.23984_wp / ener     ! [um]
      end do
      close(u)

      open(newunit=u, file=F_CPE, status='old', action='read')
      read(u, '(/)')
      do i = 1, NCPE
         read(u, *) ener, rn1, rk, e1, e2
         cpe_eV(i)   = ener
         cpe_eps1(i) = e1
         cpe_eps2(i) = e2
         cpe_wavl(i) = 1.23984_wp / ener
      end do
      close(u)

      loaded = .true.
   end subroutine load_tables


   subroutine free_diel(ener_eV, agrain, dtype, eps1, eps2)
      ! Free-electron contribution to graphite dielectric function.
      ! Form confirmed by Draine.
      real(wp),         intent(in)  :: ener_eV, agrain
      character(len=*), intent(in)  :: dtype
      real(wp),         intent(out) :: eps1, eps2
      real(wp) :: E_plasma, Tau_bulk, Veff, Tau, X, G

      if (dtype == 'CpeD03') then
         ! graphite, E perpendicular to c-axis (basal plane)
         E_plasma = 0.285_wp * sqrt(1.0_wp - 6.24e-3_wp*T_GRAIN + 3.66e-5_wp*T_GRAIN*T_GRAIN)
         Tau_bulk = 4.2e-11_wp / (1.0_wp + 0.322_wp*T_GRAIN + 1.30e-3_wp*T_GRAIN*T_GRAIN)
         Veff     = 4.5e11_wp * sqrt(1.0_wp + T_GRAIN/255.0_wp)
      else
         ! graphite, E parallel to c-axis ('CpaD03'); Tau_bulk = 3e-14 s
         ! per Draine (vs. 1.4e-14 in Draine & Lee 1984).
         E_plasma = 0.101_wp
         Tau_bulk = 3.0e-14_wp
         Veff     = 3.7e10_wp * sqrt(1.0_wp + T_GRAIN/255.0_wp)
      end if

      Tau  = 1.0_wp / (1.0_wp/Tau_bulk + Veff/agrain)
      X    = ener_eV / E_plasma
      G    = 1.0_wp / (1.518e15_wp * E_plasma * Tau)
      eps1 = -1.0_wp / (X*X + G*G)
      eps2 =  G / (X*(X*X + G*G))
   end subroutine free_diel


   subroutine eps_to_nk(eps1, eps2, n, k)
      ! eps1 here is Re(eps) − 1 (Draine convention). Handles small-eps
      ! cases to avoid catastrophic cancellation in deep-IR / hard X-ray.
      real(wp), intent(in)  :: eps1, eps2
      real(wp), intent(out) :: n, k
      real(wp) :: rr

      if ((eps1*eps1 + eps2*eps2) < 1.0e-6_wp) then
         n = 0.5_wp*eps1 - 0.125_wp*(eps1*eps1 - eps2*eps2) + 1.0_wp
         k = 0.5_wp*eps2 - 0.25_wp*eps1*eps2
      else if (eps2 < 1.0e-3_wp * abs(1.0_wp + eps1)) then
         n = sqrt(1.0_wp + eps1) * (1.0_wp + 0.125_wp*(eps2/(1.0_wp + eps1))**2)
         k = sqrt(1.0_wp + eps1) * 0.5_wp * (eps2/(1.0_wp + eps1))
      else
         rr = sqrt((1.0_wp + eps1)**2 + eps2*eps2)
         n  = sqrt(0.5_wp*(rr + 1.0_wp + eps1))
         k  = sqrt(0.5_wp*(rr - 1.0_wp - eps1))
      end if
   end subroutine eps_to_nk


   subroutine build_nk(agrain)
      ! Refresh cached (n, k) tables for `agrain`.
      real(wp), intent(in) :: agrain
      integer  :: i
      real(wp) :: e1f, e2f, e1, e2

      do i = 1, NCPA
         call free_diel(cpa_eV(i), agrain, 'CpaD03', e1f, e2f)
         e1 = cpa_eps1(i) + e1f
         e2 = cpa_eps2(i) + e2f
         call eps_to_nk(e1, e2, cpa_n(i), cpa_k(i))
      end do
      do i = 1, NCPE
         call free_diel(cpe_eV(i), agrain, 'CpeD03', e1f, e2f)
         e1 = cpe_eps1(i) + e1f
         e2 = cpe_eps2(i) + e2f
         call eps_to_nk(e1, e2, cpe_n(i), cpe_k(i))
      end do
      cached_a = agrain
   end subroutine build_nk


   subroutine q_graphite_abs(agrain, lambda, Qabs)
      ! Random-orientation-averaged Q_abs for a graphite sphere.
      ! agrain, lambda: μm. Qabs: dimensionless (C_abs / π a²).
      real(wp), intent(in)  :: agrain
      real(wp), intent(in)  :: lambda
      real(wp), intent(out) :: Qabs
      real(wp) :: x, n_pa, k_pa, n_pe, k_pe
      real(wp) :: Qext1, Qsca1, alb1, gsca1
      real(wp) :: Qabs_pa, Qabs_pe

      if (.not. loaded) call load_tables()
      if (agrain /= cached_a) call build_nk(agrain)

      ! cpa_wavl / cpe_wavl are descending (eV is ascending in file);
      ! interp1 handles both monotonic directions.
      call interp(cpa_wavl, cpa_n, lambda, n_pa)
      call interp(cpa_wavl, cpa_k, lambda, k_pa)
      call interp(cpe_wavl, cpe_n, lambda, n_pe)
      call interp(cpe_wavl, cpe_k, lambda, k_pe)

      x = 2.0_wp * PI_LOC * agrain / lambda
      call mie(n_pa, k_pa, x, Qext1, Qsca1, Qabs_pa, alb1, gsca1)
      call mie(n_pe, k_pe, x, Qext1, Qsca1, Qabs_pe, alb1, gsca1)

      Qabs = Qabs_pa/3.0_wp + 2.0_wp*Qabs_pe/3.0_wp
   end subroutine q_graphite_abs


   subroutine q_graphite_full(agrain, lambda, Qext, Qsca, Qabs, gsca)
      ! Random-orientation-averaged full Mie output for a graphite sphere.
      ! Orientation average is (1/3) E||c + (2/3) E_|_c applied to each
      ! efficiency; g is averaged scattering-weighted (the physically correct
      ! combination). agrain, lambda: um.
      real(wp), intent(in)  :: agrain, lambda
      real(wp), intent(out) :: Qext, Qsca, Qabs, gsca
      real(wp) :: x, n_pa, k_pa, n_pe, k_pe
      real(wp) :: Qext_pa, Qsca_pa, Qabs_pa, alb_pa, g_pa
      real(wp) :: Qext_pe, Qsca_pe, Qabs_pe, alb_pe, g_pe

      if (.not. loaded) call load_tables()
      if (agrain /= cached_a) call build_nk(agrain)

      call interp(cpa_wavl, cpa_n, lambda, n_pa)
      call interp(cpa_wavl, cpa_k, lambda, k_pa)
      call interp(cpe_wavl, cpe_n, lambda, n_pe)
      call interp(cpe_wavl, cpe_k, lambda, k_pe)

      x = 2.0_wp * PI_LOC * agrain / lambda
      call mie(n_pa, k_pa, x, Qext_pa, Qsca_pa, Qabs_pa, alb_pa, g_pa)
      call mie(n_pe, k_pe, x, Qext_pe, Qsca_pe, Qabs_pe, alb_pe, g_pe)

      Qabs = Qabs_pa/3.0_wp + 2.0_wp*Qabs_pe/3.0_wp
      Qsca = Qsca_pa/3.0_wp + 2.0_wp*Qsca_pe/3.0_wp
      Qext = Qabs + Qsca
      if (Qsca > 0.0_wp) then
         gsca = (Qsca_pa*g_pa/3.0_wp + 2.0_wp*Qsca_pe*g_pe/3.0_wp) / Qsca
      else
         gsca = 0.0_wp
      end if
   end subroutine q_graphite_full

end module q_graphite_mod
