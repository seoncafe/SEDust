module read_index
   ! Reader for the DH21a astrodust dielectric-function tables
   ! `index_DH21Ad_Pporo_fFe_ba`.
   !
   ! File format (described in
   ! astrodust_Draine_Hensley/astrodust_DielectricFunction/Readme.txt):
   !   line 1     : ICOMP / parameter description
   !   line 2     : column header
   !   lines 3-6336 (6334 rows):  E[eV]  Re(m)-1  Im(m)  Re(eps)-1  Im(eps)
   !
   ! Rows are in ascending E, which means descending lambda. We reverse on
   ! load so the internal arrays are in ascending lambda — that keeps
   ! bisection straightforward and matches the convention used elsewhere
   ! (DH21_wave is also ascending lambda).
   !
   ! `interp_m(lambda, nr, ki)` returns Re(m) and Im(m) at lambda [microns]
   ! by linear interpolation in log(lambda). m itself is interpolated
   ! linearly (not log) because Re(m)-1 changes sign at high energy where
   ! the material is X-ray transparent.
   !
   ! No EMA / mixing is applied here; the file already encodes the
   ! Bruggeman effective dielectric function for the requested (P, fFe, b/a).

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   implicit none
   private
   public :: load_index, interp_m, ndata, lambda, n_r, k_i

   integer, parameter :: wp = real64

   integer  :: ndata = 0
   real(wp), allocatable :: lambda(:)   ! [microns], ascending
   real(wp), allocatable :: n_r(:)      ! Re(m)   = (Re(m)-1) + 1
   real(wp), allocatable :: k_i(:)      ! Im(m)
   real(wp), allocatable :: log_lam(:)  ! log10(lambda), cached for interp

contains

   subroutine load_index(filename)
      character(len=*), intent(in) :: filename
      integer  :: u, ios, i, j
      real(wp) :: e_ev, rem1, imm, reeps1, imeps
      real(wp), allocatable :: e_buf(:), nr_buf(:), ki_buf(:)
      character(len=512) :: line

      ! First pass to count data lines (skip 2 header lines).
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') 'load_index: cannot open ', trim(filename)
         stop 1
      end if
      read(u,'(a)') line   ! header 1
      read(u,'(a)') line   ! header 2
      ndata = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         ndata = ndata + 1
      end do
      close(u)

      if (ndata < 2) then
         write(error_unit,'(a,i0,a)') 'load_index: only ', ndata, &
              ' data rows found; expected >= 6000.'
         stop 1
      end if

      allocate(e_buf(ndata), nr_buf(ndata), ki_buf(ndata))

      ! Second pass to actually read the values, in original (ascending E) order.
      open(newunit=u, file=filename, status='old', action='read')
      read(u,'(a)') line
      read(u,'(a)') line
      do i = 1, ndata
         read(u,*) e_ev, rem1, imm, reeps1, imeps
         e_buf(i)  = e_ev
         nr_buf(i) = rem1 + 1.0_wp
         ki_buf(i) = imm
      end do
      close(u)

      ! Reverse to ascending lambda. lambda[um] = 1.2398 / E[eV].
      if (allocated(lambda))  deallocate(lambda)
      if (allocated(n_r))     deallocate(n_r)
      if (allocated(k_i))     deallocate(k_i)
      if (allocated(log_lam)) deallocate(log_lam)
      allocate(lambda(ndata), n_r(ndata), k_i(ndata), log_lam(ndata))
      do i = 1, ndata
         j         = ndata - i + 1
         lambda(i) = 1.2398_wp / e_buf(j)
         n_r(i)    = nr_buf(j)
         k_i(i)    = ki_buf(j)
         log_lam(i) = log10(lambda(i))
      end do

      deallocate(e_buf, nr_buf, ki_buf)

      ! Sanity: monotonicity in lambda.
      do i = 2, ndata
         if (lambda(i) <= lambda(i-1)) then
            write(error_unit,'(a,i0)') &
               'load_index: lambda not strictly increasing at i=', i
            stop 1
         end if
      end do
   end subroutine load_index


   subroutine interp_m(lam, nr_out, ki_out)
      ! Linear interpolation of (n_r, k_i) in log(lambda).
      real(wp), intent(in)  :: lam            ! [microns]
      real(wp), intent(out) :: nr_out, ki_out
      integer  :: lo, hi, mid
      real(wp) :: t, x

      if (ndata < 2) then
         write(error_unit,'(a)') 'interp_m: load_index has not been called.'
         stop 1
      end if
      if (lam < lambda(1) .or. lam > lambda(ndata)) then
         write(error_unit,'(a,es12.4,a,es12.4,a,es12.4,a)') &
            'interp_m: lambda = ', lam, ' [um] outside table range [', &
            lambda(1), ', ', lambda(ndata), '].'
         stop 1
      end if

      x  = log10(lam)
      lo = 1
      hi = ndata
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (log_lam(mid) <= x) then
            lo = mid
         else
            hi = mid
         end if
      end do
      t      = (x - log_lam(lo)) / (log_lam(hi) - log_lam(lo))
      nr_out = (1.0_wp - t)*n_r(lo) + t*n_r(hi)
      ki_out = (1.0_wp - t)*k_i(lo) + t*k_i(hi)
   end subroutine interp_m

end module read_index
