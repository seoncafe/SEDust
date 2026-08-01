module kext_table_mod
   ! Reader and wavelength interpolant for the size-integrated extinction
   ! tables calc_kext.x writes into data/ (kext_astrodust_MW*.dat,
   ! kext_dl07_MW*.dat, kext_zubko_BARE_GR_S.dat, kext_<model>.dat).
   !
   ! These files ARE the transport optics of a dust model: extinction,
   ! absorption and scattering cross section per H nucleon, the scattering
   ! albedo and the scattering asymmetry <cos>, already integrated over the
   ! size distribution of every population.  A radiative-transfer host that
   ! links libsedust.a reads them through dust_extinction rather than
   ! re-running the size integral in the middle of a transport run.
   !
   ! FILE FORMAT.  Comment lines start with '#', blank lines are ignored, and
   ! every remaining line is one wavelength:
   !
   !   lambda[um]  albedo  <cos>  C_ext/H  C_abs/H  C_sca/H  [K_abs]  [C_polext/H]
   !
   ! Only the first six fields are read.  The trailing K_abs [cm^2/g] and
   ! dichroic C_polext/H columns are written by some products and not by
   ! others (the EUV products are scalar optics only), so the number of
   ! columns varies between files and a reader must not depend on it.
   ! Wavelength must be positive and strictly increasing.
   use constants, only: wp
   implicit none
   private
   public :: load_kext_table, tabulated_extinction_on_grid
   public :: LAM_MATCH_TOL

   ! Relative tolerance at which two wavelengths count as the SAME grid node.
   ! It has to be non-zero because a table stores lambda in decimal with
   ! 6-13 significant digits while the model carries the full binary value, so
   ! a node the two share can differ in the last digits after the round trip
   ! through the file.  1e-6 is the tolerance load_q_table already uses to
   ! decide that two rows sit at one wavelength, and the tables are spaced
   ! dln(lambda) ~ 0.012 apart -- four orders of magnitude coarser -- so it can
   ! only ever snap a round-trip difference, never two distinct nodes
   ! together.  The same tolerance widens the end points of the table for the
   ! range test, so a grid that reaches exactly as far as the table is not
   ! rejected by its own last decimal digit.
   real(wp), parameter :: LAM_MATCH_TOL = 1.0e-6_wp

contains

   subroutine load_kext_table(path, n, lam, albedo, gbar, Cext, Cabs, Csca, ok)
      ! Read one extinction table.  On success ok = .true., n is the number of
      ! data rows and the six arrays are allocated to (n).  On any failure --
      ! the file is absent, it holds fewer than two data rows, a row cannot be
      ! parsed into six reals, or the wavelengths are not positive and strictly
      ! increasing -- ok = .false., n = 0 and nothing is left allocated.
      !
      ! A failure is REPORTED, never fatal: the caller (a model builder) has to
      ! be able to distinguish a host's explicitly named file, whose absence is
      ! a configuration error, from the built-in default, whose absence only
      ! means this model serves no extinction to an RT host.
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: n
      real(wp), allocatable, intent(out) :: lam(:), albedo(:), gbar(:)
      real(wp), allocatable, intent(out) :: Cext(:), Cabs(:), Csca(:)
      logical,               intent(out) :: ok

      integer            :: u, ios, i
      character(len=512) :: line
      real(wp)           :: row(6)

      n = 0;  ok = .false.

      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return

      ! Pass 1: count the data rows.
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         n = n + 1
      end do
      if (n < 2) then
         close(u);  n = 0;  return
      end if

      allocate(lam(n), albedo(n), gbar(n), Cext(n), Cabs(n), Csca(n))

      ! Pass 2: read them.
      rewind(u)
      i = 0
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
         i = i + 1
         if (i > n) exit
         read(line, *, iostat=ios) row
         if (ios /= 0) then
            close(u)
            deallocate(lam, albedo, gbar, Cext, Cabs, Csca)
            n = 0;  return
         end if
         lam(i)    = row(1)
         albedo(i) = row(2)
         gbar(i)   = row(3)
         Cext(i)   = row(4)
         Cabs(i)   = row(5)
         Csca(i)   = row(6)
      end do
      close(u)

      if (i /= n) then
         deallocate(lam, albedo, gbar, Cext, Cabs, Csca)
         n = 0;  return
      end if

      if (lam(1) <= 0.0_wp) then
         deallocate(lam, albedo, gbar, Cext, Cabs, Csca)
         n = 0;  return
      end if
      do i = 2, n
         if (lam(i) <= lam(i-1)) then
            deallocate(lam, albedo, gbar, Cext, Cabs, Csca)
            n = 0;  return
         end if
      end do

      ok = .true.
   end subroutine load_kext_table


   subroutine tabulated_extinction_on_grid(nt, lam_t, Cext_t, Cabs_t, Csca_t, gbar_t, &
                                           lam, Cext, Cabs, Csca, gbar, ok)
      ! Put a loaded table onto a model's own wavelength grid.
      !
      ! The cross sections are interpolated linearly in log(lambda)-log(C),
      ! because extinction, absorption and scattering are positive and
      ! locally power-law-like in wavelength.  Where a bracketing table value
      ! is not positive the logarithm does not exist, and that interval is
      ! interpolated linearly in the value instead, on the same log(lambda)
      ! axis.  The asymmetry <cos> is a signed, bounded direction cosine, not
      ! a power law -- it changes sign in the far infrared of some models --
      ! so it is always linear in log(lambda).
      !
      ! A grid wavelength that coincides with a table node (within
      ! LAM_MATCH_TOL) takes that node's value unchanged, so a model whose grid
      ! IS the table grid gets the table back exactly rather than through a
      ! round trip of logarithms.
      !
      ! ok = .false. if any grid wavelength falls outside the table.  There is
      ! no extrapolation: past its ends the table says nothing about the
      ! optics, and inventing a value there would silently hand the transport a
      ! frozen boundary cross section.
      ! Allocates nothing, does no I/O and touches no module state, so an RT
      ! host may call it (through dust_extinction) from its photon threads.
      integer,  intent(in)  :: nt
      real(wp), intent(in)  :: lam_t(nt), Cext_t(nt), Cabs_t(nt), Csca_t(nt), gbar_t(nt)
      real(wp), intent(in)  :: lam(:)                       ! model grid [um]
      real(wp), intent(out) :: Cext(:), Cabs(:), Csca(:)
      real(wp), optional, intent(out) :: gbar(:)
      logical,  intent(out) :: ok

      integer  :: i, j, n
      real(wp) :: x, t, l1, l2

      ok = .false.
      n = size(lam)
      if (nt < 2) return

      do i = 1, n
         x = lam(i)
         if (x <= 0.0_wp) return
         if (x < lam_t(1)  * (1.0_wp - LAM_MATCH_TOL)) return
         if (x > lam_t(nt) * (1.0_wp + LAM_MATCH_TOL)) return

         j = bracket(lam_t, nt, x)          ! lam_t(j-1) <= x <= lam_t(j)

         if (abs(x - lam_t(j-1)) <= LAM_MATCH_TOL * lam_t(j-1)) then
            Cext(i) = Cext_t(j-1);  Cabs(i) = Cabs_t(j-1);  Csca(i) = Csca_t(j-1)
            if (present(gbar)) gbar(i) = gbar_t(j-1)
            cycle
         end if
         if (abs(x - lam_t(j)) <= LAM_MATCH_TOL * lam_t(j)) then
            Cext(i) = Cext_t(j);  Cabs(i) = Cabs_t(j);  Csca(i) = Csca_t(j)
            if (present(gbar)) gbar(i) = gbar_t(j)
            cycle
         end if

         l1 = log(lam_t(j-1));  l2 = log(lam_t(j))
         t  = (log(x) - l1) / (l2 - l1)
         Cext(i) = power_law_or_linear(Cext_t(j-1), Cext_t(j), t)
         Cabs(i) = power_law_or_linear(Cabs_t(j-1), Cabs_t(j), t)
         Csca(i) = power_law_or_linear(Csca_t(j-1), Csca_t(j), t)
         if (present(gbar)) gbar(i) = (1.0_wp - t) * gbar_t(j-1) + t * gbar_t(j)
      end do

      ok = .true.
   end subroutine tabulated_extinction_on_grid


   pure function power_law_or_linear(y1, y2, t) result(y)
      ! Interpolate a cross section between two table nodes at fractional
      ! position t in log(lambda): a power law in lambda where both nodes are
      ! positive, linear in the value where one of them is not (a cross section
      ! that underflows to zero has no logarithm, and a table may carry such a
      ! point at the far infrared end of a weakly scattering model).
      real(wp), intent(in) :: y1, y2, t
      real(wp)             :: y
      if (y1 > 0.0_wp .and. y2 > 0.0_wp) then
         y = exp((1.0_wp - t) * log(y1) + t * log(y2))
      else
         y = (1.0_wp - t) * y1 + t * y2
      end if
   end function power_law_or_linear


   pure function bracket(x, n, x0) result(j)
      ! Index j with x(j-1) <= x0 <= x(j) on an ascending grid; the end
      ! intervals are returned for an x0 sitting on (or a rounding step past)
      ! either end.
      integer,  intent(in) :: n
      real(wp), intent(in) :: x(n), x0
      integer :: j, lo, hi, mid
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (x(mid) <= x0) then
            lo = mid
         else
            hi = mid
         end if
      end do
      j = hi
   end function bracket

end module kext_table_mod
