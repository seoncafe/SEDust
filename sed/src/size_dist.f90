module size_dist_mod
   ! Reader for the HD23 astrodust size distribution
   !   data/release/size_distribution.dat
   !
   ! File format (5 columns, 4 comment-`#` header lines, ~167 data rows):
   !   a_eff[um]  dn_Ad/nH  dn_PAH/nH  f_ion  f_align
   !
   ! From the PLAN (B.1) and the file's own header:
   !   "dn_i/nH = (1/nH) * (dn/da)_i * da" — i.e. the number per H is
   !   ALREADY integrated over the local size bin. Do not multiply by da.
   !
   ! The astrodust-only path uses the Astrodust column. dn_PAH/nH, f_ion,
   ! and f_align are read into separate arrays (used by the PAH and
   ! polarization paths) and ignored by the tau-Ad and SED-Ad calculations.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   implicit none
   private
   public :: load_size_dist
   public :: n_size, a_dist, dn_ad, dn_pah, f_ion, f_align

   integer, parameter :: wp = real64

   integer  :: n_size = 0
   real(wp), allocatable :: a_dist(:)   ! [microns], ascending
   real(wp), allocatable :: dn_ad(:)    ! [1/H]  (already binned)
   real(wp), allocatable :: dn_pah(:)   ! [1/H]
   real(wp), allocatable :: f_ion(:)    ! PAH ionization fraction
   real(wp), allocatable :: f_align(:)  ! alignment fraction

contains

   subroutine load_size_dist(filename, ok)
      character(len=*),  intent(in)  :: filename
      ! Optional ok: absent -> stop on error as before; present -> return
      ! .false. (leaving the module unloaded) instead of stopping.
      logical, optional, intent(out) :: ok
      integer :: u, ios, i
      character(len=512) :: line

      if (present(ok)) ok = .true.
      ! First pass: count data lines (skip leading `#` lines).
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(ok)) then
            ok = .false.;  return
         else
            write(error_unit,'(a,a)') 'load_size_dist: cannot open ', trim(filename)
            stop 1
         end if
      end if
      n_size = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0)   cycle
         if (line(1:1) == '#')      cycle
         n_size = n_size + 1
      end do
      close(u)

      if (n_size < 2) then
         if (present(ok)) then
            n_size = 0;  ok = .false.;  return
         else
            write(error_unit,'(a,i0,a)') 'load_size_dist: ', n_size, &
               ' data rows; expected > 100.'
            stop 1
         end if
      end if

      if (allocated(a_dist))  deallocate(a_dist)
      if (allocated(dn_ad))   deallocate(dn_ad)
      if (allocated(dn_pah))  deallocate(dn_pah)
      if (allocated(f_ion))   deallocate(f_ion)
      if (allocated(f_align)) deallocate(f_align)
      allocate(a_dist(n_size), dn_ad(n_size), dn_pah(n_size), &
               f_ion(n_size), f_align(n_size))

      ! Second pass: read.
      open(newunit=u, file=filename, status='old', action='read')
      i = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0)   cycle
         if (line(1:1) == '#')      cycle
         i = i + 1
         read(line,*) a_dist(i), dn_ad(i), dn_pah(i), f_ion(i), f_align(i)
      end do
      close(u)

      do i = 2, n_size
         if (a_dist(i) <= a_dist(i-1)) then
            if (present(ok)) then
               deallocate(a_dist, dn_ad, dn_pah, f_ion, f_align)
               n_size = 0;  ok = .false.;  return
            else
               write(error_unit,'(a,i0)') &
                  'load_size_dist: a_dist not strictly ascending at i=', i
               stop 1
            end if
         end if
      end do
   end subroutine load_size_dist

end module size_dist_mod
