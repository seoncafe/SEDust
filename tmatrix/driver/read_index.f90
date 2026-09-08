module read_index
   ! Refractive-index tables m(lambda) = n + i k, read from the two text
   ! layouts Draine distributes.
   !
   ! LAYOUT 1 -- "index" files, ordered by photon energy.  Used by
   ! `index_DH21Ad_Pporo_fFe_ba` (astrodust, described in
   ! astrodust_Draine_Hensley/astrodust_DielectricFunction/Readme.txt) and by
   ! `index_amcBE_ZMCB96` (the BE amorphous carbon of Zubko et al. 1996):
   !
   !     E[eV]  Re(m)-1  Im(m)  Re(eps)-1  Im(eps)
   !
   ! with rows in ascending E, i.e. descending lambda, and lambda[um] =
   ! 1.2398 / E[eV].
   !
   ! LAYOUT 2 -- "eps" files, ordered by wavelength.  Used by `eps_suvSil`
   ! (the pre-2001 "smoothed UV" astrosilicate of Draine & Lee 1984 and
   ! Laor & Draine 1993, which Weingartner & Draine 2001 adopt):
   !
   !     w[um]  Re(eps-1)  Im(eps)  Re(m-1)  Im(m)
   !
   ! with rows in descending lambda.
   !
   ! Both layouts carry m in the same two columns of the same five, so one
   ! reader handles both and only the meaning of column 1 differs.  Header
   ! lines are recognized rather than counted: a line is data when it holds no
   ! '=' sign, five reals can be read from it, and its first value is
   ! positive.  Every header line of every file above fails at least one of
   ! those (they carry '=' or start with a word), and no data line of any of
   ! them does.
   !
   ! Both loaders leave the table in ASCENDING lambda, which keeps bisection
   ! straightforward and matches the convention used elsewhere (DH21_wave is
   ! also ascending lambda).
   !
   ! `refractive_index_at(table, lambda, nr, ki)` interpolates linearly in
   ! log(lambda).  m itself is interpolated linearly (not log) because
   ! Re(m)-1 changes sign at high energy where the material is X-ray
   ! transparent.
   !
   ! No effective-medium mixing is applied here.  The astrodust index file
   ! already encodes the Bruggeman effective dielectric function for the
   ! requested (P, fFe, b/a); mixing two tables of the kind read here is the
   ! business of module `effective_medium`.
   !
   ! Two interfaces are offered.  `refractive_index_t` holds one table and is
   ! what a caller needing several materials at once uses.  `load_index` and
   ! `interp_m` act on a single table owned by this module; they are the
   ! entry points every existing driver in this directory uses, and they read
   ! layout 1.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   implicit none
   private
   public :: refractive_index_t
   public :: load_index_energy_table, load_index_wavelength_table
   public :: refractive_index_at
   public :: load_index, interp_m

   integer, parameter :: wp = real64

   !! Photon energy of a 1 micron photon, in eV: the conversion constant the
   !! Draine index files are tabulated against.
   real(wp), parameter :: E_EV_MICRON = 1.2398_wp

   type :: refractive_index_t
      integer :: ndata = 0
      real(wp), allocatable :: lambda(:)    ! [micron], strictly ascending
      real(wp), allocatable :: n_r(:)       ! Re(m)
      real(wp), allocatable :: k_i(:)       ! Im(m)
      real(wp), allocatable :: log_lam(:)   ! log10(lambda), cached for interp
   end type refractive_index_t

   !! Table addressed by load_index / interp_m.
   type(refractive_index_t) :: single_table

contains

   subroutine load_index(filename)
      !! Layout-1 loader for the module's own table.
      character(len=*), intent(in) :: filename
      call load_index_energy_table(single_table, filename)
   end subroutine load_index


   subroutine interp_m(lam, nr_out, ki_out)
      !! Interpolation on the module's own table.
      real(wp), intent(in)  :: lam            ! [microns]
      real(wp), intent(out) :: nr_out, ki_out
      call refractive_index_at(single_table, lam, nr_out, ki_out)
   end subroutine interp_m


   subroutine load_index_energy_table(tab, filename)
      !! Layout 1: E[eV] Re(m)-1 Im(m) Re(eps)-1 Im(eps), ascending E.
      type(refractive_index_t), intent(out) :: tab
      character(len=*), intent(in) :: filename
      real(wp), allocatable :: col(:,:)
      integer :: n, i, j

      call read_five_column_table(filename, n, col)
      allocate(tab%lambda(n), tab%n_r(n), tab%k_i(n), tab%log_lam(n))
      tab%ndata = n
      ! Ascending E is descending lambda, so fill back to front.
      do i = 1, n
         j = n - i + 1
         tab%lambda(i) = E_EV_MICRON / col(1,j)
         tab%n_r(i)    = col(2,j) + 1.0_wp
         tab%k_i(i)    = col(3,j)
      end do
      call finish_table(tab, filename)
   end subroutine load_index_energy_table


   subroutine load_index_wavelength_table(tab, filename)
      !! Layout 2: w[um] Re(eps-1) Im(eps) Re(m-1) Im(m), descending lambda.
      type(refractive_index_t), intent(out) :: tab
      character(len=*), intent(in) :: filename
      real(wp), allocatable :: col(:,:)
      integer :: n, i, j

      call read_five_column_table(filename, n, col)
      allocate(tab%lambda(n), tab%n_r(n), tab%k_i(n), tab%log_lam(n))
      tab%ndata = n
      if (col(1,1) > col(1,n)) then
         do i = 1, n
            j = n - i + 1
            tab%lambda(i) = col(1,j)
            tab%n_r(i)    = col(4,j) + 1.0_wp
            tab%k_i(i)    = col(5,j)
         end do
      else
         do i = 1, n
            tab%lambda(i) = col(1,i)
            tab%n_r(i)    = col(4,i) + 1.0_wp
            tab%k_i(i)    = col(5,i)
         end do
      end if
      call finish_table(tab, filename)
   end subroutine load_index_wavelength_table


   subroutine read_five_column_table(filename, n, col)
      !! Reads every data line of a five-column Draine table into col(1:5, 1:n)
      !! in file order.  A line is data when it carries no '=' sign, yields
      !! five reals to a list-directed read, and has a positive first value;
      !! everything else is a header line and is skipped.
      character(len=*), intent(in) :: filename
      integer, intent(out) :: n
      real(wp), allocatable, intent(out) :: col(:,:)
      character(len=512) :: line
      real(wp) :: v(5)
      integer  :: u, ios, pass

      n = 0
      do pass = 1, 2
         open(newunit=u, file=filename, status='old', action='read', iostat=ios)
         if (ios /= 0) then
            write(error_unit,'(a,a)') &
               'read_five_column_table: cannot open ', trim(filename)
            stop 1
         end if
         if (pass == 2) then
            allocate(col(5,n))
            n = 0
         end if
         do
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (index(line,'=') > 0) cycle
            read(line,*,iostat=ios) v
            if (ios /= 0) cycle
            if (v(1) <= 0.0_wp) cycle
            n = n + 1
            if (pass == 2) col(:,n) = v
         end do
         close(u)
         if (n < 2) then
            write(error_unit,'(a,i0,a,a)') &
               'read_five_column_table: only ', n, ' data rows found in ', &
               trim(filename)
            stop 1
         end if
      end do
   end subroutine read_five_column_table


   subroutine finish_table(tab, filename)
      !! Caches log10(lambda) and asserts strict monotonicity.
      type(refractive_index_t), intent(inout) :: tab
      character(len=*), intent(in) :: filename
      integer :: i

      do i = 1, tab%ndata
         tab%log_lam(i) = log10(tab%lambda(i))
      end do
      do i = 2, tab%ndata
         if (tab%lambda(i) <= tab%lambda(i-1)) then
            write(error_unit,'(a,i0,a,a)') &
               'finish_table: lambda not strictly increasing at i=', i, &
               ' in ', trim(filename)
            stop 1
         end if
      end do
   end subroutine finish_table


   subroutine refractive_index_at(tab, lam, nr_out, ki_out)
      !! Linear interpolation of (n_r, k_i) in log(lambda).
      type(refractive_index_t), intent(in) :: tab
      real(wp), intent(in)  :: lam            ! [microns]
      real(wp), intent(out) :: nr_out, ki_out
      integer  :: lo, hi, mid
      real(wp) :: t, x

      if (tab%ndata < 2) then
         write(error_unit,'(a)') &
            'refractive_index_at: the table has not been loaded.'
         stop 1
      end if
      if (lam < tab%lambda(1) .or. lam > tab%lambda(tab%ndata)) then
         write(error_unit,'(a,es12.4,a,es12.4,a,es12.4,a)') &
            'refractive_index_at: lambda = ', lam, ' [um] outside table range [', &
            tab%lambda(1), ', ', tab%lambda(tab%ndata), '].'
         stop 1
      end if

      x  = log10(lam)
      lo = 1
      hi = tab%ndata
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (tab%log_lam(mid) <= x) then
            lo = mid
         else
            hi = mid
         end if
      end do
      t      = (x - tab%log_lam(lo)) / (tab%log_lam(hi) - tab%log_lam(lo))
      nr_out = (1.0_wp - t)*tab%n_r(lo) + t*tab%n_r(hi)
      ki_out = (1.0_wp - t)*tab%k_i(lo) + t*tab%k_i(hi)
   end subroutine refractive_index_at

end module read_index
