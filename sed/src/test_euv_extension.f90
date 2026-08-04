program test_euv_extension
   !====================================================================
   ! Checks for the extreme-ultraviolet extension of the astrodust model
   ! wavelength grid (the optional lam_min argument of build_astrodust).
   !
   ! Usage (run from sed/, data paths are relative to ../):
   !   ./test_euv_extension.x                       ! no EUV polarized table
   !   ./test_euv_extension.x EUV_Q EUV_WAVE        ! with one
   !
   ! EUV_Q / EUV_WAVE are a Q stream and its wavelength axis as written by
   !   ./run_q_jori.x lam L1 L2 ...
   ! in tmatrix/.  Given, they are handed to build_astrodust as qpol_euv_path /
   ! qpol_euv_wave_path; omitted, the built-in default path is used, which in a
   ! fresh tree does not exist -- and that is the state the deficit branch of
   ! checks 4 and 5 verifies.
   !
   ! Anchors:
   !   1. Grid -- an extended build starts exactly at lam_min, is strictly
   !      increasing, and its prepended points are spaced no more coarsely
   !      than the T-matrix Q table is at its short-wavelength end.
   !   2. Table block untouched -- lam(n_euv+1:) of the extended model is
   !      bit-identical to lam(:) of an unextended one, and so are Cext,
   !      Cabs, Csca, gbar and the albedo there. Extending the grid must not
   !      perturb a single wavelength the table already covered.
   !   3. Join at 0.0912 um -- how far the volume-equivalent-sphere optics
   !      just shortward of the Q table's end sit from the spheroid table
   !      just longward of it, measured against the size of one ordinary
   !      table-to-table step of the same dln(lambda). This is the quality
   !      indicator of the sphere approximation, not a pass/fail assertion,
   !      so the tolerance is deliberately loose.
   !   4. Polarized build -- an extended build with the orientation-resolved
   !      table loaded completes (the wavelength-grid guard in build_Cpol
   !      accepts the shorter table block) and reproduces the unextended
   !      polarized optics over the table block.
   !   5. The EUV band -- the alignment-weighted, size-integrated
   !      C_pol,ext / C_ext of the astrodust grains at and around the seam.
   !   6. DL07 -- the same grid extension on the model whose optics are
   !      dielectric-function Mie throughout, where the extension is a grid
   !      extension and nothing else.
   !
   ! CHECKS 4 AND 5 HAVE TWO BRANCHES, chosen from the data rather than from a
   ! switch: whether the EUV block of Cpol_ext and Cbir_ext came back non-zero.
   !
   !   EUV BLOCK ZERO -- the DEFAULT state of a fresh tree.  The polarized EUV
   !   optics are generated one wavelength at a time, on demand
   !   (tmatrix/driver/run_q_jori.f90, `lam` mode), so most models carry no EUV
   !   companion table at all; build_Cpol then says so on stderr and leaves the
   !   block at zero.  What is verified is that DOCUMENTED DEFICIT and nothing
   !   more: the block is EXACTLY zero (not small, not interpolated from the
   !   optical band), the table block is untouched by the extension, and the
   !   model still reports itself as carrying polarized optics -- the deficit
   !   is confined to the EUV wavelengths.  A zero there is an omission, NOT
   !   the physical value: a b/a = 1.4 spheroid has non-zero dichroic
   !   extinction throughout this band (see check 5's reference values).
   !
   !   EUV BLOCK NON-ZERO -- an EUV companion table was supplied (built by
   !   run_q_jori.x and passed through the two optional command-line arguments
   !   below).  The full physics is then checked: the sign is opposite to the
   !   optical band, the seam is continuous, and the alignment-weighted
   !   size-integrated fractions reproduce the direct T-matrix values at the
   !   seam, at the 20.6 eV peak and at the 100 eV end.
   !
   ! Passing an EUV table whose wavelengths are sparse (a few named nodes
   ! rather than a full axis) exercises the same branch, but the log(lambda)
   ! interpolation between distant nodes is then coarse and the reference
   ! comparisons are correspondingly approximate.
   !
   ! C_ext in the polarized fractions of check 5 is the ASTRODUST extinction
   ! alone, not the model total: only the astrodust grains align, so mixing the
   ! unaligned PAH extinction into the denominator would make the fraction a
   ! statement about the PAH abundance rather than about the grain optics.
   !
   ! Each check prints PASS/FAIL with numbers; a FAIL sets a non-zero exit.
   !====================================================================
   use constants, only: wp
   ! size_integrated_extinction, not dust_extinction: every check below
   ! compares one build of a model against another build of the same model, so
   ! it has to see the optics THIS build computed. dust_extinction serves an
   ! RT host from a precomputed table, which would answer both builds with the
   ! same file and measure nothing.
   use dust_lib,  only: dust_model_t, build_astrodust, build_dl07, &
                        size_integrated_extinction, &
                        dust_has_polarized_optics, dust_nlam
   implicit none

   character(len=*), parameter :: QTAB  = &
      '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'
   ! Our own orientation-resolved run, used instead of the implicit default
   ! because it carries the 4th (forward-amplitude real-part) block: with the
   ! 3-block release table Cbir_ext is zero everywhere and check 4 could not
   ! tell "zeroed in the EUV block" from "never built at all".
   character(len=*), parameter :: QPOL  = &
      '../tmatrix/output/q_astrodust_jori_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: QWAVE = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: QAEFF = '../data/dielectric/DH21_aeff'

   integer,  parameter :: NT_IN = 100
   real(wp), parameter :: T_LO = 2.7_wp, T_HI = 5.0e3_wp
   ! DL07 reference model: WD01 MW R_V = 3.1, q_PAH = 4.6% (index 7), U = 1.
   integer,  parameter :: SD_INDEX = 7
   real(wp), parameter :: U_ISRF   = 1.0_wp
   ! 100 eV = hc / 1.23984e-2 um: the hard end of the band a photoionization
   ! radiative transfer transports.
   real(wp), parameter :: LAM_MIN = 1.24e-2_wp

   ! The table block is copied, not recomputed, so its differences are 0.
   real(wp), parameter :: TOL_BLOCK = 1.0e-14_wp
   ! Grid endpoint: lam(1) is assigned lam_min directly, so this is exact too.
   real(wp), parameter :: TOL_GRID  = 1.0e-14_wp

   ! First-principles reference values for check 5: the alignment-weighted,
   ! size-integrated C_pol,ext / C_ext,astrodust of the b/a = 1.4 DH21 spheroid,
   ! from the direct T-matrix calculation in tmatrix/driver (see
   ! euv_polarized_optics.f90 and run_q_jori.f90's `euv` mode). The dichroism
   ! is NEGATIVE through the whole EUV band -- the opposite sign to the optical
   ! band -- and passes through a maximum in magnitude near 20.6 eV.
   !
   ! WHICH TREATMENT ABOVE x = 50 THESE REFER TO. They were measured with the
   ! T-matrix taken at face value wherever it returned IERR = 0, i.e. WITHOUT
   ! the two-setting certification that run_q_jori.f90 applies when it writes a
   ! table (cross_sections_large_x). A certified table leaves every x it cannot
   ! certify at the geometric-optics limit, where Q_pol,ext = 0 exactly, so it
   ! gives a SMALLER magnitude wherever those sizes carry weight. That gap
   ! closes toward long wavelengths, because x > 60 then means only the rare
   ! largest grains: at 20.6 eV it starts at a = 0.575 um and a certified table
   ! reproduces REF_PEAK_POL to ~0.4%, whereas at 100 eV it starts at
   ! a = 0.118 um, in the heart of the size distribution. There a point anchor
   ! is not the right test at all: the certified table measures -1.06e-4, the
   ! same integral with the T-matrix taken at face value gives -1.56e-4, and
   ! continuing the measured Q_pol,ext ~ 1/x decay through the uncertified
   ! region instead of zeroing it gives -3.14e-4. The truth lies between the
   ! certified value (a strict lower bound on the magnitude, since zeroing
   ! x > 60 can only remove signal) and that 1/x continuation. So 100 eV is
   ! checked as a BRACKET, not against a number one particular treatment
   ! happened to produce.
   real(wp), parameter :: LAM_SEAM  = 0.0912_wp     ! 13.6 eV, the table's end
   real(wp), parameter :: LAM_PEAK  = 0.0602_wp     ! 20.6 eV
   real(wp), parameter :: REF_SEAM_POL = -1.257e-3_wp
   real(wp), parameter :: REF_SEAM_BIR = -7.165e-3_wp
   real(wp), parameter :: REF_PEAK_POL = -3.641e-3_wp
   ! 100 eV bracket on |Cpol_ext| / Cext: certified lower bound .. 1/x upper.
   real(wp), parameter :: REF_100EV_LO = 1.00e-4_wp
   real(wp), parameter :: REF_100EV_HI = 3.20e-4_wp
   ! Loose enough to absorb the log(lambda) interpolation onto a grid whose
   ! spacing lam_min sets at run time, tight enough that a wrong table, a wrong
   ! sign or a missing band cannot pass.
   real(wp), parameter :: TOL_REF = 0.15_wp
   ! Seam continuity: the fractional jump allowed between the last EUV point
   ! and the table's first. The two come from different files computed on
   ! different wavelength axes, so they are not required to be identical.
   real(wp), parameter :: TOL_SEAM = 0.10_wp

   type(dust_model_t) :: m_base, m_euv, m_scal_base, m_scal_euv
   type(dust_model_t) :: m_dl_base, m_dl_euv
   integer :: nfail, st, nl_b, nl_e, n_euv
   real(wp), allocatable :: Cext_b(:), Cabs_b(:), Csca_b(:), gbar_b(:), alb_b(:)
   real(wp), allocatable :: Cext_e(:), Cabs_e(:), Csca_e(:), gbar_e(:), alb_e(:)
   real(wp), allocatable :: Cpol_b(:), Cbir_b(:), Cpol_e(:), Cbir_e(:)
   ! Astrodust-only extinction, the denominator of the polarized fractions.
   real(wp), allocatable :: Cad_b(:), Cad_e(:)
   ! Optional EUV companion table, taken from the command line (see the header).
   character(len=512) :: euv_q, euv_wave
   logical :: euv_given
   ! .true. iff the EUV block of the extended polarized build came back
   ! non-zero, which is what selects the branch of checks 4 and 5. Set once,
   ! from the data, right after those builds.
   logical :: euv_block_filled

   nfail = 0

   euv_given = .false.
   euv_q = '';  euv_wave = ''
   if (command_argument_count() >= 2) then
      call get_command_argument(1, euv_q)
      call get_command_argument(2, euv_wave)
      euv_given = .true.
   else if (command_argument_count() == 1) then
      write(*,'(a)') ' usage: ./test_euv_extension.x [EUV_Q EUV_WAVE]'
      stop 2
   end if

   write(*,'(a)') '==================================================================='
   write(*,'(a)') ' test_euv_extension'
   write(*,'(a)') '   qtable = '//QTAB
   write(*,'(a,es12.5,a)') '   lam_min = ', LAM_MIN, ' um'
   if (euv_given) then
      write(*,'(a)') '   EUV polarized table = '//trim(euv_q)
      write(*,'(a)') '   its wavelength axis = '//trim(euv_wave)
   else
      write(*,'(a)') '   EUV polarized table = (none given; the built-in default '// &
                     'path is used)'
   end if
   write(*,'(a)') '==================================================================='

   ! ---- scalar builds: unextended vs extended --------------------------
   ! Scalar mode first, so checks 1-3 exercise the grid and the unpolarized
   ! optics on their own, with no polarized table in the picture.
   call build_astrodust(m_scal_base, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st, &
                        load_polarized_optics=.false.)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: scalar unextended build failed, status = ', st
      stop 2
   end if
   call build_astrodust(m_scal_euv, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st, &
                        load_polarized_optics=.false., lam_min=LAM_MIN, euv_tmatrix=.false.)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: scalar extended build failed, status = ', st
      stop 2
   end if

   nl_b  = dust_nlam(m_scal_base)
   nl_e  = dust_nlam(m_scal_euv)
   n_euv = nl_e - nl_b

   allocate(Cext_b(nl_b), Cabs_b(nl_b), Csca_b(nl_b), gbar_b(nl_b), alb_b(nl_b))
   allocate(Cext_e(nl_e), Cabs_e(nl_e), Csca_e(nl_e), gbar_e(nl_e), alb_e(nl_e))
   call size_integrated_extinction(m_scal_base, Cext_b, Cabs_b, Csca_b, &
                                   gbar=gbar_b, albedo=alb_b)
   call size_integrated_extinction(m_scal_euv,  Cext_e, Cabs_e, Csca_e, &
                                   gbar=gbar_e, albedo=alb_e)

   call check_grid(nfail)
   call check_table_block(nfail)
   call check_join(nfail)

   ! ---- polarized builds: unextended vs extended -----------------------
   call build_astrodust(m_base, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st, &
                        qpol_path=QPOL, qpol_wave_path=QWAVE, qpol_aeff_path=QAEFF)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: polarized unextended build failed, status = ', st
      stop 2
   end if
   if (euv_given) then
      call build_astrodust(m_euv, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st, &
                           qpol_path=QPOL, qpol_wave_path=QWAVE, qpol_aeff_path=QAEFF, &
                           lam_min=LAM_MIN, qpol_euv_path=trim(euv_q), &
                           qpol_euv_wave_path=trim(euv_wave), euv_tmatrix=.false.)
   else
      call build_astrodust(m_euv, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st, &
                           qpol_path=QPOL, qpol_wave_path=QWAVE, qpol_aeff_path=QAEFF, &
                           lam_min=LAM_MIN, euv_tmatrix=.false.)
   end if
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: polarized extended build failed, status = ', st
      stop 2
   end if

   allocate(Cpol_b(nl_b), Cbir_b(nl_b), Cpol_e(nl_e), Cbir_e(nl_e))
   allocate(Cad_b(nl_b), Cad_e(nl_e))
   call size_integrated_extinction(m_base, Cext_b, Cabs_b, Csca_b, gbar=gbar_b, &
                        Cpol_ext=Cpol_b, Cbir_ext=Cbir_b, albedo=alb_b)
   call size_integrated_extinction(m_euv,  Cext_e, Cabs_e, Csca_e, gbar=gbar_e, &
                        Cpol_ext=Cpol_e, Cbir_ext=Cbir_e, albedo=alb_e)
   call astrodust_extinction(m_base, Cad_b)
   call astrodust_extinction(m_euv,  Cad_e)

   ! Which branch checks 4 and 5 take is a property of the data, not of the
   ! command line: an EUV table may have been asked for and still not have
   ! reached this grid. Anything non-zero anywhere in the block means it did.
   euv_block_filled = (maxval(abs(Cpol_e(1:n_euv))) > 0.0_wp .or. &
                       maxval(abs(Cbir_e(1:n_euv))) > 0.0_wp)
   write(*,'(a)') '-------------------------------------------------------------------'
   if (euv_block_filled) then
      write(*,'(a)') ' EUV polarized block: FILLED -- checks 4 and 5 verify its physics.'
   else
      write(*,'(a)') ' EUV polarized block: ABSENT -- checks 4 and 5 verify the '// &
                     'documented deficit'
      write(*,'(a)') '   (exact zeros below the seam, table block and polarized '// &
                     'capability intact).'
   end if

   call check_polarized(nfail)
   call check_euv_dichroism(nfail)

   ! ---- DL07 builds, last: they reset the module globals ---------------
   call build_dl07(m_dl_base, QTAB, SIZED, SD_INDEX, U_ISRF, NT_IN, T_LO, T_HI, status=st)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: DL07 unextended build failed, status = ', st
      stop 2
   end if
   call build_dl07(m_dl_euv, QTAB, SIZED, SD_INDEX, U_ISRF, NT_IN, T_LO, T_HI, &
                   status=st, lam_min=LAM_MIN)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: DL07 extended build failed, status = ', st
      stop 2
   end if
   call check_dl07(nfail)

   write(*,'(a)') '-------------------------------------------------------------------'
   if (nfail == 0) then
      write(*,'(a)') ' ALL CHECKS PASSED'
   else
      write(*,'(a,i0,a)') ' ', nfail, ' CHECK(S) FAILED'
   end if
   if (nfail /= 0) stop 1

contains

   ! ---- 1. grid: floor, monotonicity, spacing --------------------------
   subroutine check_grid(nf)
      integer, intent(inout) :: nf
      integer  :: k
      real(wp) :: d1, dln_max_ext, dln_table_1, dln_min
      logical  :: ok, mono

      d1 = abs(m_scal_euv%lam(1) - LAM_MIN) / LAM_MIN
      mono = .true.
      dln_min = huge(1.0_wp)
      do k = 2, nl_e
         if (m_scal_euv%lam(k) <= m_scal_euv%lam(k-1)) mono = .false.
         dln_min = min(dln_min, log(m_scal_euv%lam(k) / m_scal_euv%lam(k-1)))
      end do
      ! Coarsest step inside the prepended block, against the table's own
      ! first step. The extension must not be the sparser of the two.
      dln_max_ext = 0.0_wp
      do k = 2, n_euv + 1
         dln_max_ext = max(dln_max_ext, log(m_scal_euv%lam(k) / m_scal_euv%lam(k-1)))
      end do
      dln_table_1 = log(m_scal_base%lam(2) / m_scal_base%lam(1))

      ok = (d1 <= TOL_GRID .and. mono .and. dln_min > 0.0_wp .and. &
            dln_max_ext <= dln_table_1 * (1.0_wp + 1.0e-12_wp) .and. n_euv > 0)
      write(*,'(a)') ' [1] extended grid: floor at lam_min, monotone, spacing'
      write(*,'(a,i0,a,i0,a,i0)') '     NLAM: unextended = ', nl_b, &
           '   extended = ', nl_e, '   prepended n_euv = ', n_euv
      write(*,'(a,es12.5,a,es9.2)') '     lam(1) = ', m_scal_euv%lam(1), &
           ' um   rel err vs lam_min = ', d1
      write(*,'(a,l1)') '     strictly increasing = ', mono
      write(*,'(a,es10.3,a,es10.3)') '     max dln(lam) in extension = ', dln_max_ext, &
           '   table first dln = ', dln_table_1
      call verdict(ok, nf)
   end subroutine check_grid

   ! ---- 2. table block identical ---------------------------------------
   subroutine check_table_block(nf)
      integer, intent(inout) :: nf
      real(wp) :: dlam, dext, dabs, dsca, dg, dalb
      logical  :: ok
      dlam = maxreldiff(m_scal_euv%lam(n_euv+1:), m_scal_base%lam)
      dext = maxreldiff(Cext_e(n_euv+1:), Cext_b)
      dabs = maxreldiff(Cabs_e(n_euv+1:), Cabs_b)
      dsca = maxreldiff(Csca_e(n_euv+1:), Csca_b)
      dg   = maxreldiff(gbar_e(n_euv+1:), gbar_b)
      dalb = maxreldiff(alb_e(n_euv+1:),  alb_b)
      ok = (dlam <= TOL_BLOCK .and. dext <= TOL_BLOCK .and. dabs <= TOL_BLOCK .and. &
            dsca <= TOL_BLOCK .and. dg <= TOL_BLOCK .and. dalb <= TOL_BLOCK)
      write(*,'(a)') ' [2] table block unchanged by the extension (scalar build)'
      write(*,'(a,es10.2)') '     max rel |dlam|    = ', dlam
      write(*,'(a,es10.2,a,es10.2)') '     max rel |dCext|   = ', dext, &
           '   |dCabs|   = ', dabs
      write(*,'(a,es10.2,a,es10.2)') '     max rel |dCsca|   = ', dsca, &
           '   |dgbar|   = ', dg
      write(*,'(a,es10.2,a,es10.2)') '     max rel |dalbedo| = ', dalb, &
           '   tol = ', TOL_BLOCK
      call verdict(ok, nf)
   end subroutine check_table_block

   ! ---- 3. join at the Q table's short-wavelength end ------------------
   subroutine check_join(nf)
      ! j = n_euv is the last sphere (Mie, dielectric function) point;
      ! j = n_euv+1 is the table's first (spheroid, T-matrix) point.  The
      ! reference is the next table-to-table step, n_euv+1 -> n_euv+2, taken
      ! over a comparable dln(lambda): whatever that step does is the grid's
      ! own variation, and the excess over it is what the sphere shape costs.
      integer, intent(inout) :: nf
      integer  :: js, jt
      real(wp) :: r_ext, r_abs, r_sca, r_alb
      real(wp) :: t_ext, t_abs, t_sca, t_alb
      real(wp) :: dln_join, dln_tab
      logical  :: ok

      js = n_euv;  jt = n_euv + 1
      dln_join = log(m_scal_euv%lam(jt)   / m_scal_euv%lam(js))
      dln_tab  = log(m_scal_euv%lam(jt+1) / m_scal_euv%lam(jt))

      r_ext = reldiff(Cext_e(js), Cext_e(jt))
      r_abs = reldiff(Cabs_e(js), Cabs_e(jt))
      r_sca = reldiff(Csca_e(js), Csca_e(jt))
      r_alb = reldiff(alb_e(js),  alb_e(jt))
      t_ext = reldiff(Cext_e(jt), Cext_e(jt+1))
      t_abs = reldiff(Cabs_e(jt), Cabs_e(jt+1))
      t_sca = reldiff(Csca_e(jt), Csca_e(jt+1))
      t_alb = reldiff(alb_e(jt),  alb_e(jt+1))

      ! Quality indicator, not an assertion of agreement: require only that the
      ! sphere optics are finite, positive and of the same order as the table's.
      ok = (Cext_e(js) > 0.0_wp .and. Cabs_e(js) > 0.0_wp .and. Csca_e(js) > 0.0_wp &
            .and. r_ext < 1.0_wp .and. r_abs < 1.0_wp .and. r_sca < 1.0_wp)
      write(*,'(a)') ' [3] join at the Q table end: sphere (below) vs table (above)'
      write(*,'(a,es12.5,a,es12.5,a)') '     lambda: sphere = ', m_scal_euv%lam(js), &
           ' um   table = ', m_scal_euv%lam(jt), ' um'
      write(*,'(a,es10.3,a,es10.3)') '     dln(lam) across join = ', dln_join, &
           '   next table step = ', dln_tab
      write(*,'(a)') '                     across join    next table step'
      write(*,'(a,es14.3,es18.3)') '     Cext        ', r_ext, t_ext
      write(*,'(a,es14.3,es18.3)') '     Cabs        ', r_abs, t_abs
      write(*,'(a,es14.3,es18.3)') '     Csca        ', r_sca, t_sca
      write(*,'(a,es14.3,es18.3)') '     albedo      ', r_alb, t_alb
      call verdict(ok, nf)
   end subroutine check_join

   ! ---- 4. the polarized build over the two blocks ---------------------
   subroutine check_polarized(nf)
      ! The table block is asserted the same way in both branches: extending
      ! the grid must not move a wavelength the table already covered. The EUV
      ! block is asserted differently -- as physics when a companion table
      ! filled it, as an exact zero when none did.
      integer, intent(inout) :: nf
      real(wp) :: mpol_euv, mbir_euv, dpol, dbir, pol_max_v
      integer  :: kpol
      logical  :: ok, ok_block, has_b, has_e, neg_euv

      has_b = dust_has_polarized_optics(m_base)
      has_e = dust_has_polarized_optics(m_euv)
      kpol = maxloc(abs(Cpol_e(1:n_euv)), 1)
      mpol_euv = abs(Cpol_e(kpol))
      mbir_euv = maxval(abs(Cbir_e(1:n_euv)))
      dpol = maxreldiff(Cpol_e(n_euv+1:), Cpol_b)
      dbir = maxreldiff(Cbir_e(n_euv+1:), Cbir_b)
      pol_max_v = maxval(Cpol_b)
      ! Both table blocks must be non-trivial, or the agreement above would be
      ! the agreement of two all-zero arrays.
      ok_block = (has_b .and. has_e .and. pol_max_v > 0.0_wp &
                  .and. dpol <= TOL_BLOCK .and. dbir <= TOL_BLOCK &
                  .and. maxval(abs(Cpol_b)) > 0.0_wp .and. maxval(abs(Cbir_b)) > 0.0_wp)

      if (euv_block_filled) then
         ! Every EUV point must carry the sign opposite to the optical band,
         ! where the dichroism is positive: a sign slip would be a jori mix-up.
         neg_euv = all(Cpol_e(1:n_euv) < 0.0_wp)
         ok = (ok_block .and. mpol_euv > 0.0_wp .and. mbir_euv > 0.0_wp .and. neg_euv)
         write(*,'(a)') ' [4] polarized build with lam_min: EUV block is physics, '// &
              'table block kept'
      else
         ! An unfilled block must be an exact zero. A small non-zero number
         ! would mean something leaked in from the optical band -- clamped
         ! interpolation, say -- and would be a wrong value rather than a
         ! missing one.
         neg_euv = .false.
         ok = (ok_block .and. mpol_euv == 0.0_wp .and. mbir_euv == 0.0_wp)
         write(*,'(a)') ' [4] polarized build with lam_min: EUV block is an exact '// &
              'zero, table block kept'
      end if

      write(*,'(a,l1,a,l1)') '     has_polarized: unextended = ', has_b, &
           '   extended = ', has_e
      if (euv_block_filled) then
         write(*,'(a,es10.2,a,es12.5,a)') '     EUV block max|Cpol_ext| = ', mpol_euv, &
              ' cm^2/H  at ', m_euv%lam(kpol), ' um'
         write(*,'(a,f8.2,a)') '       = ', 100.0_wp*mpol_euv/pol_max_v, &
              '% of the peak |Cpol_ext| of the optical band, opposite sign'
         write(*,'(a,es10.2)') '     EUV block max|Cbir_ext| = ', mbir_euv
         write(*,'(a,l1)') '     EUV block Cpol_ext < 0 everywhere = ', neg_euv
      else
         write(*,'(a,i0,a,es12.5,a,es12.5,a)') '     EUV block: ', n_euv, &
              ' points, ', m_euv%lam(1), ' to ', m_euv%lam(n_euv), ' um'
         write(*,'(a,es10.2,a,es10.2)') '     EUV block max|Cpol_ext| = ', mpol_euv, &
              '   max|Cbir_ext| = ', mbir_euv
         write(*,'(a)') '     (required to be exactly 0 -- a KNOWN DEFICIT, not the '// &
              'physical value:'
         write(*,'(a)') '      a b/a = 1.4 spheroid is dichroic throughout this band)'
      end if
      write(*,'(a,es10.2,a,es10.2)') '     table block max rel |dCpol_ext| = ', dpol, &
           '   |dCbir_ext| = ', dbir
      write(*,'(a,es10.2,a,es10.2)') '     table block max |Cpol_ext| = ', &
           maxval(abs(Cpol_b)), '   max |Cbir_ext| = ', maxval(abs(Cbir_b))
      call verdict(ok, nf)
   end subroutine check_polarized

   ! ---- 5. the EUV dichroism against the direct T-matrix calculation ---
   subroutine check_euv_dichroism(nf)
      ! The seam values and the optical-band sign reversal are read off the
      ! MAIN table, which is present either way, so they are asserted in both
      ! branches. Only the two quantities that live below the seam -- the
      ! continuity of the join and the 20.6 eV / 100 eV fractions -- depend on
      ! the companion table, and they are replaced by the exact-zero assertion
      ! when it is absent.
      integer, intent(inout) :: nf
      integer  :: jt, js, jpk, j1059, j1072, j100
      real(wp) :: f_seam, f_below, f_peak, f_100, b_seam, jump
      real(wp) :: f_1059, f_1072
      logical  :: ok, ok_seam, zero_euv

      jt = n_euv + 1                       ! table's first (shortest) wavelength
      js = n_euv                           ! last EUV point, just below it
      f_seam  = Cpol_e(jt) / Cad_e(jt)
      b_seam  = Cbir_e(jt) / Cad_e(jt)
      f_below = Cpol_e(js) / Cad_e(js)
      jump    = abs(f_below - f_seam) / max(abs(f_seam), tiny(1.0_wp))

      jpk   = nearest_lam(m_euv%lam, LAM_PEAK)
      j100  = 1                            ! lam(1) = lam_min = 100 eV
      j1059 = nearest_lam(m_euv%lam, 0.1059_wp)
      j1072 = nearest_lam(m_euv%lam, 0.1072_wp)
      f_peak = Cpol_e(jpk)  / Cad_e(jpk)
      f_100  = Cpol_e(j100) / Cad_e(j100)
      f_1059 = Cpol_e(j1059)/ Cad_e(j1059)
      f_1072 = Cpol_e(j1072)/ Cad_e(j1072)

      ! Anchors that need only the main table.
      ok_seam = (within(f_seam, REF_SEAM_POL, TOL_REF) &
                 .and. within(b_seam, REF_SEAM_BIR, TOL_REF) &
                 .and. f_1059 < 0.0_wp .and. f_1072 > 0.0_wp)

      if (euv_block_filled) then
         ok = (ok_seam .and. jump <= TOL_SEAM &
               .and. within(f_peak, REF_PEAK_POL, TOL_REF) &
               .and. f_100 < 0.0_wp &
               .and. abs(f_100) >= REF_100EV_LO .and. abs(f_100) <= REF_100EV_HI)
         write(*,'(a)') ' [5] EUV dichroism vs the direct T-matrix values'
      else
         zero_euv = (all(Cpol_e(1:n_euv) == 0.0_wp) .and. all(Cbir_e(1:n_euv) == 0.0_wp))
         ok = (ok_seam .and. zero_euv)
         write(*,'(a)') ' [5] seam anchors from the main table; EUV band left at zero'
      end if

      write(*,'(a)') '     (Cpol_ext and Cbir_ext over the ASTRODUST extinction)'
      write(*,'(a,es12.5,a,es12.4)') '     seam        lam = ', m_euv%lam(jt), &
           ' um   Cpol_ext/Cext = ', f_seam
      write(*,'(a,es12.5,a,es12.4)') '     just below  lam = ', m_euv%lam(js), &
           ' um   Cpol_ext/Cext = ', f_below
      write(*,'(a)') '                      model         reference     rel err'
      call report(' seam   Cpol ', f_seam, REF_SEAM_POL)
      call report(' seam   Cbir ', b_seam, REF_SEAM_BIR)
      if (euv_block_filled) then
         write(*,'(a,es10.3,a,es10.3)') '     seam jump = ', jump, '   tol = ', TOL_SEAM
         call report(' 20.6eV Cpol ', f_peak, REF_PEAK_POL)
         write(*,'(a,es12.4,a,es9.2,a,es9.2,a)') '  100 eV Cpol  ', f_100, &
              '   bracket [-', REF_100EV_HI, ', -', REF_100EV_LO, ']'
      else
         write(*,'(a,l1)') '     whole EUV band exactly zero (Cpol_ext and Cbir_ext) = ', &
              zero_euv
         write(*,'(a)') '     not compared with the 20.6 eV and 100 eV references:'
         write(*,'(a,es11.3,a,es9.2,a,es9.2,a)') '       they are ', REF_PEAK_POL, &
              ' and the bracket [-', REF_100EV_HI, ', -', REF_100EV_LO, ']'
         write(*,'(a)') '       -- the band is MISSING here, not zero by physics.'
      end if
      write(*,'(a,es12.5,a,es12.4)') '     sign reversal: at ', m_euv%lam(j1059), &
           ' um   Cpol_ext/Cext = ', f_1059
      write(*,'(a,es12.5,a,es12.4)') '                    at ', m_euv%lam(j1072), &
           ' um   Cpol_ext/Cext = ', f_1072
      call verdict(ok, nf)
   end subroutine check_euv_dichroism

   subroutine report(label, v, ref)
      character(len=*), intent(in) :: label
      real(wp),         intent(in) :: v, ref
      write(*,'(a,a,2es14.4,es12.2)') '    ', label, v, ref, abs(v-ref)/abs(ref)
   end subroutine report

   logical function within(v, ref, tol)
      real(wp), intent(in) :: v, ref, tol
      within = (abs(v - ref) <= tol * abs(ref))
   end function within

   integer function nearest_lam(lam, target) result(j)
      real(wp), intent(in) :: lam(:), target
      j = minloc(abs(lam - target), 1)
   end function nearest_lam

   subroutine astrodust_extinction(m, Cad)
      ! Size-integrated extinction of the ASTRODUST population alone,
      ! sum_a dn(a) * (Cabs + Csca), the denominator the polarized fractions
      ! are quoted against. Populations 2 and 3 are the unaligned PAHs.
      type(dust_model_t), intent(in)  :: m
      real(wp),           intent(out) :: Cad(:)
      integer :: ja
      Cad = 0.0_wp
      do ja = 1, size(m%pops(1)%dn)
         Cad = Cad + m%pops(1)%dn(ja) * (m%pops(1)%Cabs(:, ja) + m%pops(1)%Csca(:, ja))
      end do
   end subroutine astrodust_extinction

   ! ---- 6. DL07: grid extension only -----------------------------------
   subroutine check_dl07(nf)
      ! The DL07 silicate and carbonaceous optics are Mie on the D03
      ! dielectric functions, valid to 6.2e-5 um, so the extension changes
      ! nothing but the grid: the table block must come back untouched and the
      ! new points must carry finite, positive cross sections.
      integer, intent(inout) :: nf
      integer  :: nb, ne, nx
      real(wp) :: dlam, dext, dabs, dsca
      real(wp), allocatable :: xb(:), ab(:), sb(:), gb(:)
      real(wp), allocatable :: xe(:), ae(:), se(:), ge(:)
      logical  :: ok

      nb = dust_nlam(m_dl_base);  ne = dust_nlam(m_dl_euv);  nx = ne - nb
      allocate(xb(nb), ab(nb), sb(nb), gb(nb), xe(ne), ae(ne), se(ne), ge(ne))
      call size_integrated_extinction(m_dl_base, xb, ab, sb, gbar=gb)
      call size_integrated_extinction(m_dl_euv,  xe, ae, se, gbar=ge)

      dlam = maxreldiff(m_dl_euv%lam(nx+1:), m_dl_base%lam)
      dext = maxreldiff(xe(nx+1:), xb)
      dabs = maxreldiff(ae(nx+1:), ab)
      dsca = maxreldiff(se(nx+1:), sb)
      ok = (nx == n_euv .and. dlam <= TOL_BLOCK .and. dext <= TOL_BLOCK .and. &
            dabs <= TOL_BLOCK .and. dsca <= TOL_BLOCK .and. &
            minval(xe(1:nx)) > 0.0_wp .and. minval(ae(1:nx)) > 0.0_wp)
      write(*,'(a)') ' [6] DL07 with lam_min: grid extension only'
      write(*,'(a,i0,a,i0,a,i0)') '     NLAM: unextended = ', nb, &
           '   extended = ', ne, '   prepended = ', nx
      write(*,'(a,es10.2,a,es10.2)') '     table block max rel |dlam| = ', dlam, &
           '   |dCext| = ', dext
      write(*,'(a,es10.2,a,es10.2)') '     table block max rel |dCabs| = ', dabs, &
           '   |dCsca| = ', dsca
      write(*,'(a,es12.4,a,es12.4)') '     at lam_min: Cext = ', xe(1), &
           ' cm^2/H   albedo = ', se(1)/xe(1)
      call verdict(ok, nf)
      deallocate(xb, ab, sb, gb, xe, ae, se, ge)
   end subroutine check_dl07

   ! ---- utilities -----------------------------------------------------
   real(wp) function maxreldiff(a, b) result(d)
      ! Maximum over the grid of |a-b| / max(|a|,|b|,tiny). Exactly 0 when the
      ! two arrays are bit-identical.
      real(wp), intent(in) :: a(:), b(:)
      integer  :: k
      real(wp) :: s
      d = 0.0_wp
      do k = 1, size(a)
         s = max(abs(a(k)), abs(b(k)), tiny(1.0_wp))
         d = max(d, abs(a(k) - b(k)) / s)
      end do
   end function maxreldiff

   real(wp) function reldiff(a, b) result(d)
      real(wp), intent(in) :: a, b
      d = abs(a - b) / max(abs(a), abs(b), tiny(1.0_wp))
   end function reldiff

   subroutine verdict(ok, nf)
      logical, intent(in)    :: ok
      integer, intent(inout) :: nf
      if (ok) then
         write(*,'(a)') '     -> PASS'
      else
         write(*,'(a)') '     -> FAIL'
         nf = nf + 1
      end if
   end subroutine verdict

end program test_euv_extension
