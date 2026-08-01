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
   ! It shows the whole host flow: build a model once, take the transport
   ! optics from it, and take one cell's emission from the same model -- so
   ! that the light the host removes from a ray and the light it puts back
   ! refer to the same grains.
   !
   ! READING THE OPTICS FROM A FILE INSTEAD.  sed/calc_kext.x writes the very
   ! same dust_extinction curves to ../data/kext_*.dat in the column order of
   ! Draine's kext_albedo tables.  Read those when the host only transports
   ! (no dust emission), or when it is not Fortran and cannot link this
   ! library; link the library when the host also needs emission, because a
   ! file cannot supply the size-resolved absorption the heating solver wants.
   ! There is no library call that reads those files -- it would only hide
   ! which of the two the host is really using.
   use constants, only: wp
   use radfield,  only: J_Mathis
   use dust_lib,  only: dust_model_t, build_astrodust, &
                        dust_emission, dust_extinction, &
                        dust_nlam, dust_n_channel
   implicit none
   character(len=*), parameter :: QTAB  = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'
   ! Shortest wavelength an ionizing-band host needs; 0.0124 um = 100 eV.
   real(wp), parameter :: LAM_MIN_EUV = 0.0124_wp
   type(dust_model_t)    :: m, m_euv
   real(wp), allocatable :: J(:), total(:), chan(:,:)
   real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), gbar(:), albedo(:)
   real(wp), allocatable :: Cext2(:), Cabs2(:), Csca2(:)
   integer :: ipk, n, st

   ! --- (1) load a model once (here: astrodust) ------------------------
   ! status is optional; with it a failed build is reported instead of
   ! stopping the process.  Its codes are listed in build_astrodust's own
   ! header.  build_dl07, build_zubko and build_from_files take the place of
   ! this call for the other models and are used the same way from here on.
   call build_astrodust(m, QTAB, SIZED, 200, 2.7_wp, 5.0e3_wp, status=st)
   if (st /= 0) then
      print '(a,i0)', ' build_astrodust failed, status = ', st
      stop 1
   end if
   n = dust_nlam(m)
   allocate(J(n), total(n), chan(n, dust_n_channel(m)))
   allocate(Cext(n), Cabs(n), Csca(n), gbar(n), albedo(n))

   ! --- (2) transport optics: what a ray sees --------------------------
   ! Cross sections per H nucleon [cm^2/H]; gbar is the scattering-weighted
   ! <cos> for the phase function and albedo = Csca/Cext.  Both are optional
   ! and are named here, so they stay correct however the argument list grows.
   call dust_extinction(m, Cext, Cabs, Csca, gbar=gbar, albedo=albedo, status=st)
   if (st /= 0) then
      print '(a,i0)', ' dust_extinction failed, status = ', st
      stop 1
   end if

   ! --- (3) one cell: assemble the local field, get emission -----------
   call J_Mathis(1.585_wp, m%lam, J)
   call dust_emission(m, J, total, chan)

   ipk = maxloc(total, 1)
   print '(a)',            ' === external RT link to libsedust.a: OK ==='
   print '(a,a,a,i0,a,i0)', '   model=', trim(m%name), '  NLAM=', n, &
                            '  n_channel=', dust_n_channel(m)
   print '(a,es11.4,a,es11.4,a)', '   grid ', m%lam(1), ' - ', m%lam(n), ' um'
   print '(a,es12.5,a,f8.5,a,f8.5)', '   V band: C_ext/H=', at_lambda(Cext, 0.55_wp), &
                            '  albedo=', at_lambda(albedo, 0.55_wp), &
                            '  <cos>=', at_lambda(gbar, 0.55_wp)
   print '(a,f7.1,a,es12.5)', '   SED peak at lam=', m%lam(ipk), ' um, lamI/NH=', total(ipk)

   ! --- (4) a host that transports ionizing radiation ------------------
   ! The astrodust grid is the T-matrix Q table's, which stops at 0.0912 um
   ! (13.6 eV).  Passing lam_min carries it down to whatever the model's own
   ! dielectric function covers; everything else about the API is unchanged,
   ! and OMITTING lam_min leaves the grid -- and every number above -- exactly
   ! as it was.  A lam_min the dielectric function cannot reach is refused
   ! through status rather than served with a frozen refractive index.
   call build_astrodust(m_euv, QTAB, SIZED, 200, 2.7_wp, 5.0e3_wp, &
                        status=st, lam_min=LAM_MIN_EUV)
   if (st /= 0) then
      print '(a,i0)', ' EUV build_astrodust failed, status = ', st
      stop 1
   end if
   allocate(Cext2(dust_nlam(m_euv)), Cabs2(dust_nlam(m_euv)), Csca2(dust_nlam(m_euv)))
   call dust_extinction(m_euv, Cext2, Cabs2, Csca2)
   print '(a,i0,a,es11.4,a,es11.4,a)', '   with lam_min: NLAM=', dust_nlam(m_euv), &
                            ', grid ', m_euv%lam(1), ' - ', m_euv%lam(dust_nlam(m_euv)), ' um'
   print '(a,es12.5)', '   C_ext/H at 100 eV (0.0124 um) = ', Cext2(1)

contains

   function at_lambda(y, lam0) result(y0)
      ! Nearest grid point of m's own wavelength axis -- enough for a printout.
      real(wp), intent(in) :: y(:), lam0
      real(wp) :: y0
      y0 = y(minloc(abs(m%lam - lam0), 1))
   end function at_lambda

end program use_dustlib
