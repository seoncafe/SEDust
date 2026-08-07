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
   !
   ! WHERE THE NUMBERS COME FROM, AND WHY THE TWO HALVES DIFFER. One
   ! dust_extinction call fills its arguments from two different places:
   !
   !   Cext, Cabs, Csca, gbar  READ FROM A TABLE -- one of the
   !                           ../data/kext_*.dat products of calc_kext.x,
   !                           attached to the model at build time and
   !                           interpolated onto m%lam.
   !   Cpol_ext, Cbir_ext      COMPUTED FROM THE MODEL'S OWN ORIENTATION-
   !                           RESOLVED OPTICS ON EVERY CALL.
   !
   ! The polarized pair cannot be tabulated, because its alignment weight
   ! f_align(a) is RUNTIME state: a host resets it cell by cell through
   ! dust_set_alignment, and a table fixed at build time could not follow. The
   ! scalars have no such freedom -- alignment is a size weight that never
   ! enters the energy balance -- so their size integral is done once, by
   ! calc_kext.x, and recorded.
   !
   ! THE CONSEQUENCE THAT CAN SURPRISE A HOST: after dust_set_alignment, a
   ! second dust_extinction call returns CHANGED Cpol_ext and Cbir_ext and
   ! UNCHANGED Cext, Cabs, Csca, gbar. That is the physics, not a stale cache:
   ! alignment does not change how much light the grains remove, only how much
   ! of it they remove preferentially in one polarization. The run below
   ! demonstrates both halves of that statement.
   !
   ! DO NOT REPLACE dust_extinction WITH size_integrated_extinction. The size
   ! integral is public under its own name, with the identical argument list,
   ! and computes ALL SIX outputs from the model's optics -- which makes it
   ! look like the way to get the scalars and the polarized pair "from one
   ! computation". It is not the recommended path, and a polarized host should
   ! not call it unless it is absolutely required to:
   !   * it is SLOW -- a population x size x wavelength triple sum on every
   !     call, and a host that resets alignment cell by cell is precisely the
   !     one tempted to put it inside the cell loop;
   !   * there is NOTHING TO RECOMPUTE on the scalar side -- the table is that
   !     same integral, already done and recorded, and a dust model's transport
   !     optics do not depend on the transport;
   !   * MIXING THE TWO ROUTES DRIFTS -- the table route comes back through the
   !     precision written to the file, and interpolates wherever the model
   !     grid and the table do not share a node, so scalars taken from the
   !     integral in one place and from the table in another leave a host
   !     inconsistent with itself.
   ! dust_extinction already does the intended thing -- table for the scalars,
   ! fresh computation for Cpol_ext and Cbir_ext. size_integrated_extinction is
   ! for WRITING a table (what calc_kext.x does), or for a host that wants to
   ! generate a product on a lam_min grid of its own.
   !
   ! PER H, NOT PER GRAM. Cext and Cpol_ext are cross sections per H nucleon,
   ! and lamI_total / lamI_pol are emission per H. A host carrying a dust mass
   ! density instead divides by dust_mass_per_H(m) [g/H], a wavelength- and
   ! temperature-independent model constant to be evaluated once (see
   ! use_dustlib.f90 for a worked call). The polarization fractions printed
   ! below are ratios, so that choice cancels out of them.
   !
   ! WHICH WAVELENGTHS ARE POLARIZED. Cpol_ext, Cbir_ext and lamI_pol are built
   ! from the orientation-resolved OPTICS table, which covers 0.0912 - 39810 um
   ! -- 1129 of the 1762 wavelengths this example runs on. The grid is the
   ! SCALAR Q table's and reaches 1.0e-4 um (12398 eV), so over the 633
   ! wavelengths below 0.0912 um the polarized pair is zero BY OMISSION, not by
   ! physics, and build_Cpol says so on stderr rather than letting the zero pass
   ! for a result. A host that needs those wavelengths polarized generates the
   ! EUV companion table (tmatrix/driver/run_q_jori.f90, `euv` mode) and passes
   ! it as qpol_euv_path / qpol_euv_wave_path. The scalar Cext / Cabs / Csca are
   ! unaffected: their extinction table covers the whole grid.
   use constants, only: wp
   use radfield,  only: J_Mathis
   use dust_lib,  only: dust_model_t, build_dust, dust_emission, &
                        dust_extinction, dust_set_alignment, dust_nlam
   implicit none
   ! The scalar example grid is the EUV companion table.  The historical
   ! non-EUV table remains the default for ordinary SED/MC runs.
   ! One data directory: build_dust resolves this model's scalar optics AND
   ! its orientation-resolved table from data_dir/astrodust/.
   character(len=*), parameter :: DATA = '../data'
   real(wp), parameter :: DEG   = acos(-1.0_wp)/180.0_wp
   ! --- cell geometry: the host's numbers, not SEDust's ---
   real(wp), parameter :: GAMMA = 90.0_wp * DEG   ! field vs. line of sight; 90 deg = max
   real(wp), parameter :: PHI   = 30.0_wp * DEG   ! position angle of projected field
   real(wp), parameter :: FTURB = 1.0_wp          ! turbulent depolarization

   type(dust_model_t)    :: m
   real(wp), allocatable :: J(:), total(:), pol(:), total2(:), pol2(:)
   real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), Cpol_ext(:)
   real(wp), allocatable :: Cext2(:), Cabs2(:), Csca2(:), Cpol2(:)
   real(wp) :: geo, jQ, jU, jI, dtot, dext
   integer  :: n, i, k, iw(3)

   ! --- load a model once ---
   ! kext_path is omitted, so the scalar outputs of dust_extinction come from
   ! /kext of this model's own product, read on the whole wavelength axis.
   call build_dust(m, 'astrodust', DATA, 200, 2.7_wp, 5.0e3_wp, include_euv=.true.)
   n = dust_nlam(m)
   allocate(J(n), total(n), pol(n), total2(n), pol2(n))
   allocate(Cext(n), Cabs(n), Csca(n), Cpol_ext(n))
   allocate(Cext2(n), Cabs2(n), Csca2(n), Cpol2(n))

   ! --- one cell: local field -> emission, polarized part included ---
   call J_Mathis(1.585_wp, m%lam, J)
   call dust_emission(m, J, total, lamI_pol=pol)

   ! --- same model object -> opacity on the same grid, dichroic part included ---
   ! Cext/Cabs/Csca are served from the table; Cpol_ext is computed here and now
   ! with the alignment the model currently holds.
   call dust_extinction(m, Cext, Cabs, Csca, Cpol_ext=Cpol_ext)

   iw = [ilam_near(0.55_wp), ilam_near(154.0_wp), ilam_near(850.0_wp)]

   print '(a)', ' === SEDust polarization interface ==='
   print '(a,a,a,f6.3,a,f7.4,a,f5.2)', '   model=', trim(m%name), &
         '   alignment: f_max=', m%align_fmax, ' a_align=', m%align_a, &
         ' um alpha=', m%align_alpha
   print '(a,f5.1,a,f5.1,a,f4.2)', '   geometry: gamma=', GAMMA/DEG, &
         ' deg  phi=', PHI/DEG, ' deg  F_turb=', FTURB

   ! --- the host's job: project onto the sky ---
   ! p_ext below divides a freshly computed Cpol_ext by a tabulated Cext. Both
   ! describe the same grains and the same size distribution, so the ratio is
   ! the model's dichroic fraction; only the route differs.
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

   ! The extinction side of the same statement, and the asymmetry a host has to
   ! expect: this second call recomputes Cpol_ext with the new f_align, while
   ! Cext/Cabs/Csca are re-read from the same table rows and come back unchanged.
   call dust_extinction(m, Cext2, Cabs2, Csca2, Cpol_ext=Cpol2)
   dext = maxval(abs(Cext2 - Cext))

   print '(a)', ' --- f_max halved, no re-solve ---'
   print '(a,es10.3,a,es10.3)', '   max|lamI_total change|=', dtot, &
         '   max lamI_total=', maxval(total)
   print '(a,es10.3,a,es10.3)', '   max|C_ext change|     =', dext, &
         '   max C_ext      =', maxval(Cext)
   do k = 1, 3
      i = iw(k)
      print '(a,f9.2,a,f10.6,a,f10.6)', '   lam=', m%lam(i), &
            ' um   lamI_pol ratio=', pol2(i)/pol(i), &
            '   Cpol_ext ratio=', Cpol2(i)/Cpol_ext(i)
   end do

   print '(a)', ' --- division of labor ---'
   print '(a)', '   SEDust: size integral, f_align(a) weight, lamI_pol, Cpol_ext'
   print '(a)', '   RT code: sin^2(gamma), F_turb, position angle, Stokes transport'
   print '(a)', '   scalars from /kext of the model product; Cpol_ext, Cbir_ext computed per call'

contains

   ! index of the grid wavelength closest to lam0 [um]
   integer function ilam_near(lam0)
      real(wp), intent(in) :: lam0
      ilam_near = minloc(abs(m%lam - lam0), 1)
   end function ilam_near

end program use_dustlib_pol
