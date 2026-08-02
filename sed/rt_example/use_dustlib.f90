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
   ! WHERE THE TRANSPORT OPTICS COME FROM.  sed/calc_kext.x calls
   ! size_integrated_extinction to integrate the model's populations over grain
   ! size, and writes the resulting curves to ../data/kext_*.dat in the column
   ! order of Draine's kext_albedo tables.  dust_extinction serves one of those
   ! files back on the model's own wavelength grid: the builder loads it (see
   ! kext_path at step 1), because the paths are relative to sed/ and this
   ! program -- like a real host -- may change directory once the model is
   ! built.  A host that is not Fortran, or that only transports and never
   ! emits, can read the file directly instead; link the library when the host
   ! also needs emission, because a file cannot supply the size-resolved
   ! absorption the heating solver wants.
   !
   ! size_integrated_extinction, the routine that does that integral, is public
   ! as well and takes the same arguments as dust_extinction.  An RT host
   ! should NOT call it unless it is absolutely required to; dust_extinction is
   ! the normal path.  Three reasons:
   !   1. It is slow.  Every call re-does the triple sum over populations,
   !      grain sizes and wavelengths.  Nothing about it belongs in a cell loop.
   !   2. There is nothing to gain.  The table is that same integral, performed
   !      once and written down, and a model's transport optics do not depend on
   !      the transport -- so recomputing them returns the number that is
   !      already on disk.
   !   3. Mixing the two paths within one host can make its own numbers
   !      disagree.  The table path comes back through the precision the file
   !      was written with, and through interpolation wherever the model grid
   !      misses a table node.  Taking some wavelengths from the table and
   !      others from a fresh integral breaks that consistency.
   ! The legitimate use is MAKING a table -- what calc_kext.x does, and what a
   ! host would do to produce one matching a lam_min grid of its own.  Anywhere
   ! else, call dust_extinction.  This example never calls the integral.
   !
   ! PER H NUCLEON OR PER GRAM.  The cross sections dust_extinction returns are
   ! per H nucleon.  A host whose cells carry a dust mass density instead of an
   ! H column density divides them by dust_mass_per_H(m) [g/H] to get a mass
   ! opacity [cm^2/g]; the V band line this program prints is exactly that
   ! division.  The denominator is the same constant the K_abs column of every
   ! ../data/kext_*.dat is normalized by, so a host reading the file and a host
   ! linking the library arrive at the same opacity.
   use constants, only: wp
   use radfield,  only: J_Mathis
   use dust_lib,  only: dust_model_t, build_astrodust, &
                        dust_emission, dust_extinction, dust_mass_per_H, &
                        dust_nlam, dust_n_channel
   implicit none
   character(len=*), parameter :: QTAB  = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'
   ! Size-integrated extinction table this host transports with; see step (1).
   character(len=*), parameter :: KEXT_MW = '../data/kext_astrodust_MW.dat'
   ! Shortest wavelength an ionizing-band host needs; 0.0124 um = 100 eV.
   real(wp), parameter :: LAM_MIN_EUV = 0.0124_wp
   type(dust_model_t)    :: m, m_euv
   real(wp), allocatable :: J(:), total(:), chan(:,:)
   real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), gbar(:), albedo(:)
   real(wp), allocatable :: Cext2(:), Cabs2(:), Csca2(:)
   real(wp) :: Mdust_H
   integer :: ipk, n, st

   ! --- (1) load a model once (here: astrodust) ------------------------
   ! status is optional; with it a failed build is reported instead of
   ! stopping the process.  Its codes are listed in build_astrodust's own
   ! header.  build_dl07, build_zubko and build_from_files take the place of
   ! this call for the other models and are used the same way from here on.
   !
   ! kext_path names the extinction table this model will serve, and naming one
   ! is a contract: a file that cannot be read fails the build (status 5 for
   ! this builder), which is what a host wants when the file it was configured
   ! with is missing.  The file named here is the size integral computed on the
   ! T-matrix Q table's own grid, which is exactly the range this model covers
   ! -- the choice for a host that stops at the Lyman limit.  Omitting the
   ! argument falls to the model's default, ../data/kext_astrodust_MW_euv.dat:
   ! the same size integral carried further into the ultraviolet, so it CONTAINS
   ! this range and serves the plain grid just as well, which is why step (4)
   ! leaves it to the default.  A default that cannot be read is only an offer
   ! -- it leaves the model with no extinction to serve, and dust_extinction
   ! then reports status 2 -- but it does not fail the build.
   !
   ! One practical difference between the two files, since this program names
   ! the narrow one: kext_astrodust_MW.dat is the frozen regression product and
   ! carries 8 significant digits, while the EUV table carries 13.  The optics
   ! are the same integral either way -- the V band numbers below are identical
   ! to the printed precision -- but a host that wants every digit the library
   ! can give it should take the default.
   call build_astrodust(m, QTAB, SIZED, 200, 2.7_wp, 5.0e3_wp, status=st, &
                        kext_path=KEXT_MW)
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
   ! The curve comes from the table build_astrodust loaded above; status 2 would
   ! mean no table was found, and status 3 that m%lam runs outside it.
   call dust_extinction(m, Cext, Cabs, Csca, gbar=gbar, albedo=albedo, status=st)
   if (st /= 0) then
      print '(a,i0)', ' dust_extinction failed, status = ', st
      stop 1
   end if

   ! Those cross sections are per H nucleon.  A host whose cells carry a dust
   ! mass density rather than an H column density needs the mass opacity
   !    kappa_abs [cm^2/g] = (C_abs/H) / (M_dust/N_H)
   ! instead, and dust_mass_per_H supplies the denominator from the model's own
   ! size distribution and grain densities.  It depends on neither wavelength
   ! nor temperature nor radiation field, so it is evaluated ONCE here, not per
   ! cell and not per wavelength.
   Mdust_H = dust_mass_per_H(m)

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
   print '(a,es11.4,a)', '   M_dust/N_H=', Mdust_H, ' g/H  (one constant per model)'
   print '(a,es12.5,a,es12.5,a)', '   V band: C_abs/H=', at_lambda(Cabs, 0.55_wp), &
                            ' cm^2/H  ->  kappa_abs=', at_lambda(Cabs, 0.55_wp)/Mdust_H, &
                            ' cm^2/g'
   print '(a,f7.1,a,es12.5)', '   SED peak at lam=', m%lam(ipk), ' um, lamI/NH=', total(ipk)

   ! --- (4) a host that transports ionizing radiation ------------------
   ! The astrodust grid is the T-matrix Q table's, which stops at 0.0912 um
   ! (13.6 eV).  Passing lam_min carries it down to whatever the model's own
   ! dielectric function covers; everything else about the API is unchanged,
   ! and OMITTING lam_min leaves the grid -- and every number above -- exactly
   ! as it was.  A lam_min the dielectric function cannot reach is refused
   ! through status rather than served with a frozen refractive index.
   !
   ! No kext_path this time: the default IS the table that reaches down here.
   ! Handing this model the KEXT_MW of step (1) would build without complaint
   ! -- the file reads fine -- and then dust_extinction would refuse the call
   ! with status 3, because the grid now runs below anything that table
   ! tabulates and the library extrapolates no optics.
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
