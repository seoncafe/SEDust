program check_data_root
   ! Can a host build every model from a working directory that is not a
   ! SEDust sed/, naming its data with an ABSOLUTE path -- and does a bad path
   ! come back as a status rather than as a runtime abort?
   !
   !   ./check_data_root.x <absolute data dir>
   !
   ! Both halves used to fail. data_dir resolved the model directories but not
   ! the dielectric functions, which were compile-time paths relative to the
   ! working directory, so a host that moved its data separated a model from
   ! the optical constants its optics were computed on.  And those opens
   ! carried no iostat, so the Fortran runtime killed the image before the
   ! builder could set a code -- which an MPI host cannot recover from, one
   ! rank dying uncollected tearing down the job before any rank can report.
   !
   ! Run it from anywhere EXCEPT sed/, or it proves nothing: from sed/ the old
   ! relative paths resolve by accident.
   !
   !   PASS  every model builds with status 0 from the absolute directory, and
   !         a directory that does not exist returns non-zero from every one
   !   FAIL  any model returns non-zero on the good directory
   !   (an abort here is a failure too, and a louder one)
   use constants, only: wp
   use dust_lib,  only: dust_model_t, build_dust, dust_extinction
   implicit none

   character(len=512) :: ddir
   character(len=*), parameter :: MODELS(4) = &
        [character(len=9) :: 'astrodust', 'dl07', 'mrn', 'zubko']
   type(dust_model_t) :: m
   character(len=160) :: msg
   real(wp), allocatable :: Cext(:), Cabs(:), Csca(:)
   integer :: i, st, kst, nbad

   if (command_argument_count() < 1) then
      write(*,'(a)') ' usage: check_data_root.x <absolute data dir>'
      stop 1
   end if
   call get_command_argument(1, ddir)
   if (ddir(1:1) /= '/') then
      write(*,'(a)') ' check_data_root: give an ABSOLUTE data directory'
      stop 1
   end if

   nbad = 0
   write(*,'(a,a)') ' data_dir  = ', trim(ddir)
   write(*,'(a)')   ' ---- every model, from this working directory ----'
   do i = 1, size(MODELS)
      call build_dust(m, trim(MODELS(i)), trim(ddir), 40, 2.7_wp, 3.0e3_wp, &
                      include_euv = .false., status = st, message = msg)
      if (st /= 0) then
         nbad = nbad + 1
         write(*,'(a,a12,a,i0,a,a)') '   ', trim(MODELS(i)), '  status = ', st, &
                                     '   ', trim(msg)
         cycle
      end if
      allocate(Cext(m%NLAM), Cabs(m%NLAM), Csca(m%NLAM))
      call dust_extinction(m, Cext, Cabs, Csca, status = kst)
      write(*,'(a,a12,a,i0,a,es12.5,a,es12.5,a,i0)') '   ', trim(MODELS(i)), &
           '  status = ', st, '   NLAM = ', real(m%NLAM, wp), &
           '   lambda(1) = ', m%lam(1), ' um   kext status = ', kst
      deallocate(Cext, Cabs, Csca)
   end do

   write(*,'(a)') ' ---- a directory that is not there: status, not an abort ----'
   do i = 1, size(MODELS)
      call build_dust(m, trim(MODELS(i)), '/nonexistent/sedust/data', 40, &
                      2.7_wp, 3.0e3_wp, include_euv = .false., status = st, &
                      message = msg)
      write(*,'(a,a12,a,i0,a,a)') '   ', trim(MODELS(i)), '  status = ', st, &
                                  '   ', trim(msg)
      if (st == 0) nbad = nbad + 1
   end do

   if (nbad == 0) then
      write(*,'(a)') ' PASS'
   else
      write(*,'(a,i0,a)') ' FAIL (', nbad, ' case(s))'
      stop 1
   end if
end program check_data_root
