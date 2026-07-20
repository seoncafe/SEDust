module q_silicate_mod
   ! Q_abs(a, lambda) for amorphous astronomical silicate spheres.
   !
   ! Draine (2003) "astrosilicate" optical constants (index_silD03),
   ! fed through Mie theory. Silicate is isotropic (single dielectric,
   ! no orientation average) and has no free-electron contribution, so
   ! the refractive index depends only on wavelength, not grain radius
   ! (unlike graphite). Mirrors the structure of q_graphite_mod but
   ! reads (n, k) directly from the index table.
   !
   ! Reads the Draine (2003) silicate dielectric function and drives the
   ! silicate Mie call.

   use constants, only: wp
   use sed_mathlib,   only: interp
   use mie_mod,   only: mie
   implicit none
   private
   public :: q_silicate_abs
   public :: q_silicate_full

   character(len=*), parameter :: F_SIL = '../data/dielectric/index_silD03'
   integer,  parameter :: NSIL = 837
   real(wp), parameter :: PI_LOC = 3.141592653589793238462643383279502884197_wp

   logical  :: loaded = .false.
   real(wp) :: sil_eV(NSIL), sil_n(NSIL), sil_k(NSIL), sil_wavl(NSIL)

contains

   subroutine load_tables()
      ! index_silD03 columns: E[eV]  Re(n)-1  Im(n)  Re(eps)-1  Im(eps)
      ! (2 header lines). Silicate n,k read directly: n = 1 + col2, k = col3.
      integer  :: i, u
      real(wp) :: ener, rn1, rk, e1, e2

      open(newunit=u, file=F_SIL, status='old', action='read')
      read(u, '(/)')
      do i = 1, NSIL
         read(u, *) ener, rn1, rk, e1, e2
         sil_eV(i)   = ener
         sil_n(i)    = 1.0_wp + rn1
         sil_k(i)    = rk
         sil_wavl(i) = 1.23984_wp / ener     ! [um]
      end do
      close(u)
      loaded = .true.
   end subroutine load_tables


   subroutine q_silicate_abs(agrain, lambda, Qabs)
      ! Q_abs for a silicate sphere. agrain, lambda: um. Qabs: C_abs/(pi a^2).
      real(wp), intent(in)  :: agrain, lambda
      real(wp), intent(out) :: Qabs
      real(wp) :: x, nr, ki, Qext1, Qsca1, alb1, gsca1

      if (.not. loaded) call load_tables()

      ! sil_wavl is descending (eV ascending in file); interp handles both.
      call interp(sil_wavl, sil_n, lambda, nr)
      call interp(sil_wavl, sil_k, lambda, ki)

      x = 2.0_wp * PI_LOC * agrain / lambda
      call mie(nr, ki, x, Qext1, Qsca1, Qabs, alb1, gsca1)
   end subroutine q_silicate_abs


   subroutine q_silicate_full(agrain, lambda, Qext, Qsca, Qabs, gsca)
      ! Full Mie output (extinction, scattering, absorption efficiencies and
      ! scattering asymmetry g) for a silicate sphere. agrain, lambda: um.
      ! Same dielectric path as q_silicate_abs -- just keeps every Mie return.
      real(wp), intent(in)  :: agrain, lambda
      real(wp), intent(out) :: Qext, Qsca, Qabs, gsca
      real(wp) :: x, nr, ki, alb1

      if (.not. loaded) call load_tables()

      call interp(sil_wavl, sil_n, lambda, nr)
      call interp(sil_wavl, sil_k, lambda, ki)

      x = 2.0_wp * PI_LOC * agrain / lambda
      call mie(nr, ki, x, Qext, Qsca, Qabs, alb1, gsca)
   end subroutine q_silicate_full

end module q_silicate_mod
