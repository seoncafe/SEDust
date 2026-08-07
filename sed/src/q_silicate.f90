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
   use sed_paths, only: sed_data_path
   implicit none
   private
   public :: q_silicate_abs
   public :: q_silicate_full
   public :: silicate_index_lambda_range
   public :: silicate_index_available

   character(len=*), parameter :: F_SIL = 'dielectric/index_silD03'
   integer,  parameter :: NSIL = 837
   real(wp), parameter :: PI_LOC = 3.141592653589793238462643383279502884197_wp

   logical  :: loaded = .false.
   ! Whether the load succeeded.  Separate from `loaded`, which only records
   ! that it was attempted: a caller must be able to learn that the file was
   ! not there and refuse the build, instead of the open aborting the process.
   logical  :: load_ok = .false.
   real(wp) :: sil_eV(NSIL), sil_n(NSIL), sil_k(NSIL), sil_wavl(NSIL)

contains

   subroutine silicate_index_lambda_range(lam_lo, lam_hi)
      ! Shortest and longest wavelength [um] the D03 astrosilicate dielectric
      ! function covers.  Outside it `interp` returns the boundary value, i.e.
      ! a CONSTANT (n, k), so a caller whose wavelength grid runs past this
      ! range must refuse rather than let the frozen index pass as physics.
      real(wp), intent(out) :: lam_lo, lam_hi

      if (.not. loaded) call load_tables()
      lam_lo = minval(sil_wavl)
      lam_hi = maxval(sil_wavl)
   end subroutine silicate_index_lambda_range


   logical function silicate_index_available()
      ! Is the D03 astrosilicate dielectric function readable where the data
      ! root points?  Triggers the load, so a builder calls this once before
      ! committing to a model and reports through its own status.
      if (.not. loaded) call load_tables()
      silicate_index_available = load_ok
   end function silicate_index_available


   subroutine load_tables()
      ! index_silD03 columns: E[eV]  Re(n)-1  Im(n)  Re(eps)-1  Im(eps)
      ! (2 header lines). Silicate n,k read directly: n = 1 + col2, k = col3.
      integer  :: i, u, ios
      real(wp) :: ener, rn1, rk, e1, e2

      loaded  = .true.
      load_ok = .false.
      sil_eV = 0.0_wp;  sil_n = 1.0_wp;  sil_k = 0.0_wp;  sil_wavl = 1.0_wp

      open(newunit=u, file=sed_data_path(F_SIL), status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read(u, '(/)', iostat=ios)
      if (ios /= 0) then;  close(u);  return;  end if
      do i = 1, NSIL
         read(u, *, iostat=ios) ener, rn1, rk, e1, e2
         if (ios /= 0) then;  close(u);  return;  end if
         sil_eV(i)   = ener
         sil_n(i)    = 1.0_wp + rn1
         sil_k(i)    = rk
         sil_wavl(i) = 1.23984_wp / ener     ! [um]
      end do
      close(u)
      load_ok = .true.
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
