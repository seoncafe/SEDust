program use_dustlib_pol
   ! Polarization example for an EXTERNAL Fortran code (a 3D polarized RT
   ! driver) linking the dust-emission library. It is NOT part of the sed
   ! build -- it is compiled separately against libsedust.a + the .mod search
   ! path, exactly as an RT code would.  Build and run it from sed/, because
   ! the data paths below are relative to that directory:
   !
   !   cd sed
   !   make libsedust.a
   !   gfortran -I. rt_example/use_dustlib_pol.f90 libsedust.a -fopenmp -o use_dustlib_pol.x
   !   ./use_dustlib_pol.x
   !
   ! What it shows: where SEDust stops and the RT code starts. SEDust returns
   ! the INTRINSIC polarized emission lamI_pol and dichroic extinction
   ! Cpol_ext -- the size integral and the alignment weight f_align(a) are
   ! already in them. The line-of-sight geometry is the host's:
   !
   !   j_Q = lamI_pol * sin^2(gamma) * F_turb * cos(2*phi)
   !   j_U = lamI_pol * sin^2(gamma) * F_turb * sin(2*phi)
   !   j_I = lamI_total                                  (geometry-independent)
   !
   ! gamma is the angle between the LOCAL FIELD and the LINE OF SIGHT (not a
   ! sky-plane angle), phi the position angle of the projected field, F_turb a
   ! turbulent depolarization factor. The same sin^2(gamma)*F_turb multiplies
   ! Cpol_ext when the extinction matrix is assembled.
   use constants, only: wp
   use radfield,  only: J_Mathis
   use dust_lib,  only: dust_model_t, build_astrodust, dust_emission, &
                        dust_extinction, dust_set_alignment, dust_nlam
   implicit none
   character(len=*), parameter :: QTAB  = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'
   real(wp), parameter :: DEG   = acos(-1.0_wp)/180.0_wp
   ! --- cell geometry: the host's numbers, not SEDust's ---
   real(wp), parameter :: GAMMA = 90.0_wp * DEG   ! field vs. line of sight; 90 deg = max
   real(wp), parameter :: PHI   = 30.0_wp * DEG   ! position angle of projected field
   real(wp), parameter :: FTURB = 1.0_wp          ! turbulent depolarization

   type(dust_model_t)    :: m
   real(wp), allocatable :: J(:), total(:), pol(:), total2(:), pol2(:)
   real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), Cpol_ext(:)
   real(wp) :: geo, jQ, jU, jI, dtot
   integer  :: n, i, k, iw(3)

   ! --- load a model once ---
   call build_astrodust(m, QTAB, SIZED, 200, 2.7_wp, 5.0e3_wp)
   n = dust_nlam(m)
   allocate(J(n), total(n), pol(n), total2(n), pol2(n))
   allocate(Cext(n), Cabs(n), Csca(n), Cpol_ext(n))

   ! --- one cell: local field -> emission, polarized part included ---
   call J_Mathis(1.585_wp, m%lam, J)
   call dust_emission(m, J, total, lamI_pol=pol)

   ! --- same model object -> opacity on the same grid, dichroic part included ---
   call dust_extinction(m, Cext, Cabs, Csca, Cpol_ext=Cpol_ext)

   iw = [ilam_near(0.55_wp), ilam_near(154.0_wp), ilam_near(850.0_wp)]

   print '(a)', ' === SEDust polarization interface ==='
   print '(a,a,a,f6.3,a,f7.4,a,f5.2)', '   model=', trim(m%name), &
         '   alignment: f_max=', m%align_fmax, ' a_align=', m%align_a, &
         ' um alpha=', m%align_alpha
   print '(a,f5.1,a,f5.1,a,f4.2)', '   geometry: gamma=', GAMMA/DEG, &
         ' deg  phi=', PHI/DEG, ' deg  F_turb=', FTURB

   ! --- the host's job: project onto the sky ---
   geo = sin(GAMMA)**2 * FTURB
   print '(a)', '   lam[um]     lamI_total      j_Q          j_U        p_emis    p_ext'
   do k = 1, 3
      i  = iw(k)
      jI = total(i)
      jQ = pol(i) * geo * cos(2.0_wp*PHI)
      jU = pol(i) * geo * sin(2.0_wp*PHI)
      print '(f10.2,4es13.4,f9.4)', m%lam(i), jI, jQ, jU, &
            sqrt(jQ**2 + jU**2)/jI, geo*Cpol_ext(i)/Cext(i)
   end do

   ! --- alignment is a size weight outside the temperature solve: halving
   !     f_max halves the polarization and leaves lamI_total untouched ---
   call dust_set_alignment(m, 0.5_wp*m%align_fmax, m%align_a, m%align_alpha)
   call dust_emission(m, J, total2, lamI_pol=pol2)
   dtot = maxval(abs(total2 - total))

   print '(a)', ' --- f_max halved, no re-solve ---'
   print '(a,es10.3,a,es10.3)', '   max|lamI_total change|=', dtot, &
         '   max lamI_total=', maxval(total)
   do k = 1, 3
      i = iw(k)
      print '(a,f9.2,a,f10.6)', '   lam=', m%lam(i), ' um   pol ratio=', pol2(i)/pol(i)
   end do

   print '(a)', ' --- division of labor ---'
   print '(a)', '   SEDust: size integral, f_align(a) weight, lamI_pol, Cpol_ext'
   print '(a)', '   RT code: sin^2(gamma), F_turb, position angle, Stokes transport'

contains

   ! index of the grid wavelength closest to lam0 [um]
   integer function ilam_near(lam0)
      real(wp), intent(in) :: lam0
      ilam_near = minloc(abs(m%lam - lam0), 1)
   end function ilam_near

end program use_dustlib_pol
