program use_dustlib
   ! Minimal example of an EXTERNAL Fortran code (a 3D RT driver) linking the
   ! dust-emission library. It is NOT part of the sed build -- it is compiled
   ! separately against libsedust.a and its .mod search path, exactly as an RT
   ! code would.  Build and run it from sed/, because the data paths below are
   ! relative to that directory:
   !
   !   cd sed
   !   make libsedust.a
   !   gfortran -I. rt_example/use_dustlib.f90 -L. -lsedust -fopenmp \
   !            -o use_dustlib.x
   !   ./use_dustlib.x
   !
   ! ./build_lib.sh is the same library with its archive and .mod files in
   ! sed/lib/ -- the form an RT tree carries next to its own copy of sed/src/;
   ! link that one with -Ilib -Llib -lsedust.  Either way it is one archive and
   ! one -I, and neither of them contains the T-matrix.
   !
   ! WHAT THE T-MATRIX BUILD IS FOR, AND WHY THIS PROGRAM DOES NOT USE IT.
   ! WITH_TMATRIX=1, on the Makefile or on the script, puts the T-matrix and
   ! the oblate-spheroid optics of the astrodust extreme-ultraviolet band into
   ! that same archive, and a host that wants that band solved on the spheroid
   ! then calls
   !       use euv_astrodust_tmatrix, only: use_tmatrix_euv_band_optics
   !       call use_tmatrix_euv_band_optics()
   ! once before build_astrodust.  This example deliberately uses the shipped
   ! EUV companion table, so it does not ask for an on-the-fly band.  The
   ! astrodust EUV Q table covers 1.0e-4 - 3.981e4 um, i.e. from 12398 eV out past
   ! the sub-millimetre, so build_astrodust computes no wavelength outside it
   ! whatever lam_min it is given -- step (4) -- and the DL07 ionizing band is
   ! Mie on the D03 dielectric functions rather than a T-matrix -- step (5).
   ! The route is kept for the case it was written for: a Q table that stops
   ! shortward of what the host transports, or a model whose dielectric data
   ! reaches past its own table.  This program links the plain archive.
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
   use dust_lib,  only: dust_model_t, build_astrodust, build_dl07, &
                        dust_emission, dust_extinction, dust_mass_per_H, &
                        dust_nlam, dust_n_channel
   implicit none
   ! This example is the ionizing-band host path; ordinary SED/MC drivers use
   ! the historical non-EUV file without the `_euv` suffix.
   character(len=*), parameter :: QTAB  = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400_euv.dat'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'
   ! Size-integrated extinction table this host transports with; see step (1).
   character(len=*), parameter :: KEXT_MW = '../data/kext_astrodust_MW_euv.dat'
   ! Thermal-table grid: it fixes H(T,a) and kappB, which the extinction never
   ! touches, but both builders need one.
   integer,  parameter :: NT_IN = 200
   real(wp), parameter :: T_LO = 2.7_wp, T_HI = 5.0e3_wp
   ! A wavelength in the ionizing band to report at; 0.0124 um = 100 eV.
   real(wp), parameter :: LAM_100EV = 0.0124_wp
   ! DL07 at step (5): the WD01 Milky-Way R_V = 3.1, b_C = 6e-5 distribution
   ! (Draine's "60" model) in the MMP83 field.
   integer,  parameter :: SD_INDEX_DL07 = 7
   real(wp), parameter :: U_ISRF_DL07   = 1.0_wp
   ! Grid floor asked of that model.  The D03 dielectric functions stop at
   ! 6.1992e-5 um (20 keV); asking for exactly that value puts the request on
   ! the rounding boundary of the refusal test, so stand a little above it.
   real(wp), parameter :: LAM_MIN_DL07 = 6.21e-5_wp
   type(dust_model_t)    :: m, m_dl07
   real(wp), allocatable :: J(:), total(:), chan(:,:)
   real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), gbar(:), albedo(:)
   real(wp), allocatable :: Cext2(:), Cabs2(:), Csca2(:)
   real(wp) :: Mdust_H
   integer :: ipk, n, n2, j100, st

   ! --- (1) load a model once (here: astrodust) ------------------------
   ! status is optional; with it a failed build is reported instead of
   ! stopping the process.  Its codes are listed in build_astrodust's own
   ! header.  build_dl07, build_zubko and build_from_files take the place of
   ! this call for the other models and are used the same way from here on.
   !
   ! kext_path names the extinction table this model will serve, and naming one
   ! is a contract: a file that cannot be read fails the build (status 5 for
   ! this builder), which is what a host wants when the file it was configured
   ! with is missing.  Omitting the argument falls to the model's default,
   ! ../data/kext_astrodust_MW_euv.dat.  A default that cannot be read is only
   ! an offer -- it leaves the model with no extinction to serve, and
   ! dust_extinction then reports status 2 -- but it does not fail the build.
   !
   ! The named table is the EUV product written from the same 1762-point Q
   ! table.  The historical 1129-point product remains available for an
   ! ordinary non-ionizing host; this example intentionally exercises the
   ! wider ionizing grid.
   call build_astrodust(m, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st, &
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
   ! Emission comes from the ACTIVE model -- the one most recently built --
   ! because the solver reads the module-global working set the builder set up.
   ! So every dust_emission call for this model belongs here, before step (5)
   ! builds another one.
   call J_Mathis(1.585_wp, m%lam, J)
   call dust_emission(m, J, total, chan)

   ipk = maxloc(total, 1)
   print '(a)',            ' === external RT link to libsedust.a: OK ==='
   print '(a,a,a,i0,a,i0)', '   model=', trim(m%name), '  NLAM=', n, &
                            '  n_channel=', dust_n_channel(m)
   print '(a,es11.4,a,es11.4,a)', '   grid ', m%lam(1), ' - ', m%lam(n), ' um'
   print '(a,es12.5,a,f8.5,a,f8.5)', '   V band: C_ext/H=', at_lambda(m%lam, Cext, 0.55_wp), &
                            '  albedo=', at_lambda(m%lam, albedo, 0.55_wp), &
                            '  <cos>=', at_lambda(m%lam, gbar, 0.55_wp)
   print '(a,es11.4,a)', '   M_dust/N_H=', Mdust_H, ' g/H  (one constant per model)'
   print '(a,es12.5,a,es12.5,a)', '   V band: C_abs/H=', at_lambda(m%lam, Cabs, 0.55_wp), &
                            ' cm^2/H  ->  kappa_abs=', at_lambda(m%lam, Cabs, 0.55_wp)/Mdust_H, &
                            ' cm^2/g'
   print '(a,f7.1,a,es12.5)', '   SED peak at lam=', m%lam(ipk), ' um, lamI/NH=', total(ipk)

   ! --- (4) a host that transports ionizing radiation ------------------
   ! For astrodust there is nothing to do: steps (1) to (3) already cover that
   ! band.  The model's wavelength grid is the T-matrix Q table's, and that
   ! table runs to 1.0e-4 um (12398 eV); the extinction table named at step (1)
   ! covers the same range, so dust_extinction serves it there too.  The line
   ! below is read straight out of the step (2) arrays -- no second model, no
   ! extra argument.  It prints the grid wavelength it used, which is the
   ! nearest node to 100 eV rather than 100 eV itself.
   !
   ! lam_min, the argument that carries a model's grid below its Q table, can
   ! therefore do nothing here.  There is no value that would: at or above
   ! 1.0e-4 um the request already lies inside the table and no wavelength is
   ! prepended, and below it the build is refused, because the DH21 astrodust
   ! dielectric function the prepended band would have to be computed from
   ! itself stops at 1.000032e-4 um.  (Which refusal depends on the build: with
   ! the plain archive this program links it is status 6, no spheroid optics
   ! registered for the EUV band; with euv_tmatrix = .false. it is status 4,
   ! lam_min below the dielectric function.)  Step (5) is where lam_min still
   ! does the work it was written for.
   j100 = nearest_index(m%lam, LAM_100EV)
   print '(a,es11.4,a,es12.5,a)', '   ionizing band, no lam_min needed: C_ext/H at ', &
                            m%lam(j100), ' um =', Cext(j100), ' cm^2/H'

   ! --- (5) lam_min, on the model whose grid it still extends -----------
   ! DL07's optics are Mie on the Draine (2003) dielectric functions at every
   ! wavelength, and those reach 6.1992e-5 um (20 keV) -- past the Q table,
   ! which supplies this model's wavelength grid and nothing else.  A lam_min
   ! shorter than the table's 1.0e-4 um therefore does prepend wavelengths
   ! here, and the optics there are Mie on the dielectric function rather than
   ! an extrapolation of the table.  A lam_min those functions cannot reach is
   ! refused through status 4 rather than served with a frozen refractive
   ! index.  Everything else about the API is unchanged, and omitting lam_min
   ! leaves the grid on the table's own 1762 wavelengths.
   !
   ! No kext_path this time: the default, ../data/kext_dl07_MW_euv.dat, is the
   ! size integral computed on exactly this extended grid.  Handing this model
   ! ../data/kext_dl07_MW.dat instead would build without complaint -- the file
   ! reads fine -- and then dust_extinction would refuse the call with status
   ! 3, because that table starts at the Q table's 1.0e-4 um while the grid now
   ! runs below it, and the library extrapolates no optics past a table's end.
   ! Widening a model's grid means naming a table that covers the wider grid.
   call build_dl07(m_dl07, QTAB, SIZED, SD_INDEX_DL07, U_ISRF_DL07, &
                   NT_IN, T_LO, T_HI, status=st, lam_min=LAM_MIN_DL07)
   if (st /= 0) then
      print '(a,i0)', ' build_dl07 failed, status = ', st
      stop 1
   end if
   n2 = dust_nlam(m_dl07)
   allocate(Cext2(n2), Cabs2(n2), Csca2(n2))
   call dust_extinction(m_dl07, Cext2, Cabs2, Csca2, status=st)
   if (st /= 0) then
      print '(a,i0)', ' dl07 dust_extinction failed, status = ', st
      stop 1
   end if
   print '(a,a,a,i0,a,es11.4,a,es11.4,a)', '   model=', trim(m_dl07%name), &
                            ' with lam_min: NLAM=', n2, &
                            ', grid ', m_dl07%lam(1), ' - ', m_dl07%lam(n2), ' um'
   print '(a,es11.4,a,es12.5,a)', '   C_ext/H at the requested floor ', m_dl07%lam(1), &
                            ' um =', Cext2(1), ' cm^2/H'

contains

   integer function nearest_index(lam, lam0) result(j)
      ! Index of the grid point nearest lam0 -- enough for a printout, and the
      ! caller prints lam(j) so that no line can name a wavelength the number
      ! does not belong to.
      real(wp), intent(in) :: lam(:), lam0
      j = minloc(abs(lam - lam0), 1)
   end function nearest_index

   function at_lambda(lam, y, lam0) result(y0)
      real(wp), intent(in) :: lam(:), y(:), lam0
      real(wp) :: y0
      y0 = y(nearest_index(lam, lam0))
   end function at_lambda

end program use_dustlib
