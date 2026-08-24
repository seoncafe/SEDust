program test_size_distribution
   ! Compare hd23_size_distribution against the HD23 release table
   ! data/release/size_distribution.dat (a, dn_Ad/nH, dn_PAH/nH, f_ion,
   ! f_align; dn already integrated over each bin).  Passes when every
   ! column agrees to the table's own rounding: 2e-3 on any entry, 1e-4 on
   ! the size-integrated area and volume of astrodust.
   use, intrinsic :: iso_fortran_env, only: real64
   use size_dist_mod
   implicit none
   integer, parameter :: wp = real64
   character(len=*), parameter :: F_SD = '../data/release/size_distribution.dat'
   real(wp), allocatable :: t(:,:)
   real(wp) :: mx(5), sa2(2), sa3(2)
   character(len=512) :: line
   integer :: u, ios, n, i
   logical :: fail

   open(newunit=u, file=F_SD, status='old', action='read', iostat=ios)
   if (ios /= 0) then
      write(*,'(a)') ' test_size_distribution: cannot open '//F_SD;  stop 1
   end if
   n = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. adjustl(line(1:1)) == '#') cycle
      n = n + 1
   end do
   rewind(u)
   allocate(t(5, n))
   i = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0 .or. adjustl(line(1:1)) == '#') cycle
      i = i + 1
      read(line,*) t(:, i)
   end do
   close(u)

   call hd23_size_distribution()
   fail = n_size /= n
   if (fail) write(*,'(a,i0,a,i0)') ' size count: routine ', n_size, ', table ', n

   mx = 0;  sa2 = 0;  sa3 = 0
   do i = 1, min(n, n_size)
      mx(1) = max(mx(1), abs(a_dist(i)/t(1,i) - 1))
      if (t(2,i) > 0 .and. dn_ad(i) > 0)  mx(2) = max(mx(2), abs(dn_ad(i)/t(2,i) - 1))
      if (t(3,i) > 0 .and. dn_pah(i) > 0) mx(3) = max(mx(3), abs(dn_pah(i)/t(3,i) - 1))
      mx(4) = max(mx(4), abs(f_ion(i)/t(4,i) - 1))
      mx(5) = max(mx(5), abs(f_align(i)/t(5,i) - 1))
      sa2 = sa2 + [t(1,i)**2*t(2,i), a_dist(i)**2*dn_ad(i)]
      sa3 = sa3 + [t(1,i)**3*t(2,i), a_dist(i)**3*dn_ad(i)]
   end do
   write(*,'(a)')        ' max |routine/table - 1|'
   write(*,'(a,es10.2)') '   a        ', mx(1)
   write(*,'(a,es10.2)') '   dn_Ad    ', mx(2)
   write(*,'(a,es10.2)') '   dn_PAH   ', mx(3)
   write(*,'(a,es10.2)') '   f_ion    ', mx(4)
   write(*,'(a,es10.2)') '   f_align  ', mx(5)
   write(*,'(a,es10.2)') ' astrodust sum a^2 dn, routine/table - 1: ', sa2(2)/sa2(1) - 1
   write(*,'(a,es10.2)') ' astrodust sum a^3 dn, routine/table - 1: ', sa3(2)/sa3(1) - 1
   fail = fail .or. any(mx > 2e-3_wp) .or. abs(sa2(2)/sa2(1) - 1) > 1e-4_wp &
               .or. abs(sa3(2)/sa3(1) - 1) > 1e-4_wp
   if (fail) then
      write(*,'(a)') ' FAIL';  stop 1
   end if
   write(*,'(a)') ' PASS'
end program test_size_distribution
