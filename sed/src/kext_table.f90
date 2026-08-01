module kext_table_mod
   ! Reader for the size-integrated extinction curves that calc_kext.x writes
   ! under data/ (kext_astrodust_MW_euv.dat, kext_dl07_MW_euv.dat,
   ! kext_zubko_BARE_GR_S.dat, ...).  Such a file is a '#'-commented header
   ! followed by one row per wavelength, in ascending wavelength:
   !
   !   lambda[um]  albedo  <cos>  C_ext/H  C_abs/H  C_sca/H  [K_abs]
   !
   ! Cross sections are cm^2 per H nucleon.  The trailing K_abs [cm^2/g] column
   ! appears only when the model supplied its dust mass per H; it is a
   ! normalization of C_abs/H and is not read here, so both the six- and the
   ! seven-column products are accepted.
   !
   ! Nothing here stops the run.  A missing or malformed file comes back as
   ! ok = .false., so that the caller can turn it into a build status an RT host
   ! is able to recover from.
   use constants, only: wp
   implicit none
   private
   public :: load_kext_table

contains

   subroutine load_kext_table(path, n, lam, albedo, gbar, Cext, Cabs, Csca, ok)
      ! Read the whole curve.  The data rows are counted first and the arrays
      ! allocated to that count, so there is no built-in maximum length.
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: n                ! wavelengths read
      real(wp), allocatable, intent(out) :: lam(:)           ! [um], strictly ascending
      real(wp), allocatable, intent(out) :: albedo(:)        ! C_sca/C_ext as written
      real(wp), allocatable, intent(out) :: gbar(:)          ! <cos>, may be negative
      real(wp), allocatable, intent(out) :: Cext(:), Cabs(:), Csca(:)   ! [cm^2/H]
      ! .false. when the file cannot be opened, carries fewer than two data
      ! rows, has a row that does not parse into six reals, or is not strictly
      ! ascending in wavelength.  The arrays are then left unallocated and n = 0.
      logical,               intent(out) :: ok
      character(len=512) :: line
      integer :: u, ios, i

      n = 0;  ok = .false.
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) return

      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         if (carries_data(line)) n = n + 1
      end do
      if (n < 2) then
         close(u);  n = 0;  return
      end if

      rewind(u)
      allocate(lam(n), albedo(n), gbar(n), Cext(n), Cabs(n), Csca(n))
      i = 0
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         if (.not. carries_data(line)) cycle
         if (i >= n) exit
         i = i + 1
         read(line, *, iostat=ios) lam(i), albedo(i), gbar(i), Cext(i), Cabs(i), Csca(i)
         if (ios /= 0) then
            close(u);  call discard();  return
         end if
      end do
      close(u)
      if (i /= n) then
         call discard();  return
      end if

      ! Ascending order is what the interpolation onto a model grid assumes, and
      ! a repeated wavelength would make the bracket interval zero-wide.
      do i = 2, n
         if (lam(i) <= lam(i-1)) then
            call discard();  return
         end if
      end do
      ok = .true.

   contains

      subroutine discard()
         deallocate(lam, albedo, gbar, Cext, Cabs, Csca)
         n = 0
      end subroutine discard

   end subroutine load_kext_table


   pure logical function carries_data(line)
      ! A row carries data unless it is blank or a '#' comment.
      character(len=*), intent(in) :: line
      character(len=len(line)) :: s
      s = adjustl(line)
      carries_data = len_trim(s) > 0
      if (carries_data) carries_data = s(1:1) /= '#'
   end function carries_data

end module kext_table_mod
