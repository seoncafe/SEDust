program use_dustlib
   ! Minimal example of an EXTERNAL Fortran code (a 3D RT driver) linking the
   ! dust-emission library. It is NOT part of the sed build -- it is compiled
   ! separately against libsedust.a + the .mod search path, exactly as an RT
   ! code would.  Build and run it from sed/, because the data paths below are
   ! relative to that directory:
   !
   !   cd sed
   !   make libsedust.a
   !   gfortran -I. rt_example/use_dustlib.f90 libsedust.a -fopenmp -o use_dustlib.x
   !   ./use_dustlib.x
   !
   ! It loads a model once, then computes dust emission for one cell's field.
   use constants, only: wp
   use radfield,  only: J_Mathis
   use dust_lib,  only: dust_model_t, build_astrodust, build_dl07, build_zubko, &
                        dust_emission, dust_nlam, dust_n_channel
   implicit none
   character(len=*), parameter :: QTAB  = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'
   type(dust_model_t)    :: m
   real(wp), allocatable :: J(:), total(:), chan(:,:)
   integer :: ipk

   ! --- load a model once (here: astrodust) ---
   call build_astrodust(m, QTAB, SIZED, 200, 2.7_wp, 5.0e3_wp)
   allocate(J(dust_nlam(m)), total(dust_nlam(m)), chan(dust_nlam(m), dust_n_channel(m)))

   ! --- one cell: assemble the local field, get emission ---
   call J_Mathis(1.585_wp, m%lam, J)
   call dust_emission(m, J, total, chan)

   ipk = maxloc(total, 1)
   print '(a)',            ' === external RT link to libsedust.a: OK ==='
   print '(a,a,a,i0,a,i0)', '   model=', trim(m%name), '  NLAM=', dust_nlam(m), &
                            '  n_channel=', dust_n_channel(m)
   print '(a,f7.1,a,es12.5)', '   SED peak at lam=', m%lam(ipk), ' um, lamI/NH=', total(ipk)
end program use_dustlib
