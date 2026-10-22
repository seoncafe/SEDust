program test_planck
   use constants, only: wp
   use radfield, only: bbody, calc_bbody, planck_grid_t, prepare_planck_grid
   implicit none
   integer, parameter :: n=1200
   real(wp), parameter :: hc_kb_um = 6.62606957e-34_wp * 2.99792458e8_wp / &
                                     1.3806488e-23_wp * 1.0e6_wp
   real(wp) :: lam(n), reference(n), measured(n), t
   type(planck_grid_t) :: grid
   integer :: i,j, bad
   do i=1,n
      lam(i)=1.0e-4_wp*exp(real(i-1,wp)*log(1.0e8_wp)/real(n-1,wp))
   end do
   call prepare_planck_grid(lam,grid)
   bad=0
   do j=1,200
      t=2.0_wp*exp(real(j-1,wp)*log(5.0e4_wp)/199.0_wp)
      do i=1,n
         reference(i)=bbody(t,lam(i))
      end do
      call calc_bbody(t,lam,measured,grid)
      do i=1,n
         if (transfer(reference(i),0_8) /= transfer(measured(i),0_8)) then
            bad=bad+1
            if (bad<=5) print *, 'DIFFERENCE',i,j,reference(i),measured(i)
         end if
      end do
   end do
   call calc_bbody(1.0e-5_wp,lam,measured,grid)
   if (any(measured /= 0.0_wp)) bad=bad+1
   call calc_bbody(0.0_wp,lam,measured,grid)
   if (any(measured /= 0.0_wp)) bad=bad+1
   lam(1:6)=[-1.0_wp, 0.0_wp, 1.0e-4_wp, 2.0_wp, 1.0e4_wp, 1.0e-3_wp]
   call prepare_planck_grid(lam,grid)
   do j=1,2
      if (j==1) then
         t=2.725_wp
      else
         t=1.0e5_wp
      end if
      do i=1,n
         reference(i)=bbody(t,lam(i))
      end do
      call calc_bbody(t,lam,measured,grid)
      do i=1,n
         if (transfer(reference(i),0_8) /= transfer(measured(i),0_8)) bad=bad+1
      end do
   end do
   lam(1:4)=[hc_kb_um/(700.001_wp*2.725_wp), &
             hc_kb_um/(699.999_wp*2.725_wp), &
             hc_kb_um/(1.001e-4_wp*1.0e5_wp), &
             hc_kb_um/(0.999e-4_wp*1.0e5_wp)]
   call prepare_planck_grid(lam,grid)
   do j=1,2
      if (j==1) then
         t=2.725_wp
      else
         t=1.0e5_wp
      end if
      do i=1,n
         reference(i)=bbody(t,lam(i))
      end do
      call calc_bbody(t,lam,measured,grid)
      do i=1,n
         if (transfer(reference(i),0_8) /= transfer(measured(i),0_8)) bad=bad+1
      end do
   end do
   print *, 'BIT_DIFFERENCES=',bad
   if (bad/=0) stop 1
end program test_planck
