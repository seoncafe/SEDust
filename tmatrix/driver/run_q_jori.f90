program run_q_jori
   ! Precomputes the orientation-resolved optics (Q_ext, Q_abs, Q_sca for
   ! jori = 1, 2, 3) of the DH21 astrodust spheroid on its native
   ! (a_eff, lambda) grid, and writes them in the exact format that
   ! sed/src/q_table_jori.f90::load_q_table_jori reads.
   !
   ! Each (a_eff, lambda) node is handed to oriented_cross_sections
   ! (tmatrix/driver/tmatrix_oriented.f90), the shared first-principles core
   ! that selects the size-parameter regime and returns the three orientations.
   !
   ! THREE WAYS TO SET THE WAVELENGTHS.  Without the leading `euv` keyword the
   ! sweep runs on the DH21 axis ../data/dielectric/DH21_wave (1129 points,
   ! 0.0912-39810 um), producing the table the SED model's optical/infrared band
   ! reads.  With `euv` it runs on ../data/dielectric/DH21_wave_euv, the
   ! log-spaced extension SHORTWARD of 0.0912 um (13.6 eV) whose last node is
   ! the seam itself.  With `lam` / `lamfile` the wavelengths are NAMED
   ! directly, one number per value in microns, and only those are computed.
   ! Everything else -- regime thresholds, material, tolerances, file format --
   ! is identical, so the products differ only in the wavelengths they cover.
   !
   ! WHY A NAMED-WAVELENGTH MODE.  Polarized radiative transfer is run at a
   ! handful of wavelengths, not on a whole axis, and the scattering-matrix
   ! products already work that way (run_scatmat_aligned.x takes wavelengths on
   ! its command line).  `range` cannot serve that purpose: its JW1/JW2 are
   ! INDICES into a fixed axis, so it can only reach wavelengths that axis
   ! already carries.  `lam` takes the physical wavelength instead, computes the
   ! orientation-resolved optics there for all 169 radii, and writes both the Q
   ! stream and the companion wavelength-axis file, which is the pair
   ! sed/src/q_table_jori.f90::load_q_table_jori_euv reads.
   !
   ! LARGE SIZE PARAMETERS.  In `euv`, `lam` and `lamfile` the T-matrix is also
   ! attempted above x = 50, under the two-setting certification described at
   ! cross_sections_large_x; a node that fails it falls to geometric optics.
   ! The plain and `range` sweeps keep the older x > 50 -> geometric optics rule
   ! so that they reproduce the shipped table.
   !
   ! WHY THE EUV AXIS STOPS AT 0.0124 um (100 eV).  Above x = 60 the only
   ! description left is geometric_optics_limit, whose absorption is the
   ! opaque-grain Fresnel surface integral: it needs the chord optical depth
   ! 4 Im(m) x to be large.  For DH21 astrodust that holds through 100 eV
   ! (4 Im(m) x = 3.7 at x = 50, i.e. 2.6% transmission) and fails below it --
   ! 1.7 at 150 eV, 0.85 at 200 eV, 0.03 at 1 keV, where the grain is X-ray
   ! transparent and neither the extinction paradox nor the opaque Fresnel
   ! absorption applies.  Carrying the axis further down would fill the table
   ! with a limit known to be invalid there, so the floor is set where the
   ! physics of every branch used is still valid.
   !
   ! Usage:
   !   ./run_q_jori.x                    ! full sweep, 169 x 1129 points
   !   ./run_q_jori.x test               ! ~7 x 7 sample + full-sweep time estimate
   !   ./run_q_jori.x range JW1 JW2      ! full a range, jw in [JW1, JW2]
   !   ./run_q_jori.x merge FILE ...     ! assemble range outputs into the full file
   !   ./run_q_jori.x euv [ ... ]        ! any of the above on the EUV axis
   !   ./run_q_jori.x lam L1 [L2 ...] [ja=JA1:JA2] [tag=NAME]
   !   ./run_q_jori.x lamfile PATH   [ja=JA1:JA2] [tag=NAME]
   !   ./run_q_jori.x lammerge STEM FILE [FILE ...]
   !
   ! `lam` takes the wavelengths [um] as arguments; `lamfile` reads them from
   ! PATH, one per line ('#' comments and blank lines ignored, several values
   ! per line allowed).  They are sorted ascending and must be distinct.  Both
   ! write the Q stream AND the companion wavelength-axis file that
   ! load_q_table_jori_euv needs, under a stem that names the wavelengths (or
   ! the explicit tag=NAME).
   !
   ! ja=JA1:JA2 restricts one invocation to the radii JA1..JA2 of the 169-point
   ! DH21_aeff axis, so several PROCESSES can share one wavelength; `lammerge`
   ! then assembles their files.  Threads cannot: Mishchenko's solver hands the
   ! converged T-matrix to AMPL through COMMON /TMAT/ and keeps its working
   ! storage in further COMMON blocks (~39 MB), two of them blank, so the core
   ! is not re-entrant.  Separate processes each get their own copy.  A worked
   ! example on 72 cores, 3 radii per process:
   !
   !   for k in $(seq 0 56); do
   !     lo=$((k*3+1)); hi=$((lo+2)); [ $hi -gt 169 ] && hi=169
   !     ./run_q_jori.x lam 0.0602 ja=$lo:$hi &
   !   done; wait
   !   ./run_q_jori.x lammerge output/q_astrodust_jori_P0.20_Fe0.00_1.400.lam6.020E-02 \
   !        output/q_astrodust_jori_P0.20_Fe0.00_1.400.lam6.020E-02.ja*.dat
   !
   ! The split changes nothing but the wall time: every (lambda, a_eff) node is
   ! solved independently, and a merged file has been checked byte for byte
   ! against a single-process one.
   !
   ! CHOOSE THE WINDOWS SMALL.  The cost is not spread over the radii: a
   ! T-matrix solve grows like NMAX^4 with NMAX ~ x, so the few radii that sit
   ! just below the x = 60 ceiling carry most of the work while everything at
   ! x < 1 or x > 60 is closed-form and instant.  Contiguous 3-radius windows
   ! measured 34m48s of processor time in 9m22s of wall time for the three
   ! wavelengths 0.0124, 0.0602 and 0.0912 um -- a speedup of 3.7, not 57,
   ! because one window held the expensive radii.  One radius per window
   ! (ja=$k:$k, 169 processes) lets the operating system interleave them, and
   ! the wall time then falls to the cost of the single most expensive node.
   !
   ! Output (text, ASCII), in tmatrix/output/:
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.dat              (full / merge)
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.test.dat         (test, diagnostic columns)
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.jwJW1-JW2.dat    (range)
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.TAG.dat          (lam / lamfile / lammerge)
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.TAG.wave         (its wavelength axis)
   !   q_astrodust_jori_P0.20_Fe0.00_1.400.TAG.jaJA1-JA2.dat   (lam with ja=)
   ! and, in EUV mode, the same names with `_euv` after `jori`.  The default TAG
   ! is `lamL` for a single wavelength and `lamLMIN_LMAX_nN` for several.
   !
   ! A .TAG.dat + .TAG.wave pair is read by
   !   call load_q_table_jori_euv(q_file, wave_file, '../data/dielectric/DH21_aeff')
   ! and reaches the SED model through sed_init's qpol_euv_path /
   ! qpol_euv_wave_path.  That reader interpolates in log(lambda), so it needs at
   ! least two wavelengths; a one-wavelength file is still a valid product for
   ! reading directly, but it cannot fill a band.
   !
   ! Stream format of the full and range files (see the 12-line header written
   ! below and q_table_jori.f90:198-223): 12 header lines, then free-format
   !   do iq = 1, 4          ! iq = 1 Q_ext, 2 Q_abs, 3 Q_sca, 4 Q_re
   !     do jori = 1, 3
   !       do jw = 1, NW
   !         one record of NA = 169 a_eff values
   ! The 4th block Q_re is the real part of the forward-scattering amplitude
   ! response (birefringence twin of Q_ext); q_table_jori.f90 forms the
   ! birefringence 0.5*(Q_re(jori=3)-Q_re(jori=2)) from it.  A file written by
   ! an older code without this block (3 blocks only) still loads there.
   ! A range file writes only the jw in [JW1, JW2] records but keeps the same
   ! (iq, jori, jw) nesting.  Because jw is the innermost of the three loops,
   ! jw-window files are NOT concatenable into the full stream with a plain
   ! cat: the full stream needs all NW jw records contiguous inside each
   ! (iq, jori) block, whereas cat would interleave the blocks.  Use the
   ! `merge` mode, which reads the range files and re-emits them in full
   ! (iq, jori, jw) order.
   !
   ! A ja-window file (lam with ja=JA1:JA2) keeps every jw record but shortens
   ! each one to the a_eff columns JA1..JA2, so it is a column slice rather than
   ! a record slice; `lammerge` puts the slices back side by side.

   use, intrinsic :: iso_fortran_env, only: real64, int64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use constants, only: wp
   use read_index, only: load_index, interp_m
   use tmatrix_oriented, only: oriented_cross_sections, tmatrix_oriented_cross
   use asymptotic_optics, only: geometric_optics_limit
   implicit none

   ! Reference parameters (HD23 best fit).  Paths are relative to tmatrix/,
   ! where the Makefile drops run_q_jori.x.
   character(len=*), parameter :: f_aeff  = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: f_wave  = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_stem     = &
      'output/q_astrodust_jori_P0.20_Fe0.00_1.400'
   character(len=*), parameter :: f_out_full = f_stem//'.dat'
   character(len=*), parameter :: f_out_test = f_stem//'.test.dat'
   ! EUV band: the log-spaced axis shortward of the DH21 axis's 0.0912 um first
   ! node, whose own last node is that seam.
   character(len=*), parameter :: f_wave_euv = '../data/dielectric/DH21_wave_euv'
   character(len=*), parameter :: f_stem_euv = &
      'output/q_astrodust_jori_euv_P0.20_Fe0.00_1.400'
   character(len=*), parameter :: f_out_full_euv = f_stem_euv//'.dat'
   character(len=*), parameter :: f_out_test_euv = f_stem_euv//'.test.dat'

   integer, parameter :: NA = 169, NORI = 3
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2
   integer,  parameter :: NP_OBL  = -1            ! oblate spheroid, Mishchenko convention
   real(wp), parameter :: X_SMALL = 0.1_wp
   real(wp), parameter :: X_LARGE = 50.0_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   ! ---- large-x certification (euv, lam and lamfile) -----------------------
   ! Largest x at which a T-matrix solve is attempted beyond X_LARGE.  The
   ! Mishchenko storage limit NPN1 = 100 caps the multipole order and produces
   ! a hard stop near x ~ 83, but the certification below already rejects every
   ! node above x ~ 57 (measured at the seam: x = 51.7 and 54.7 pass, 57.9,
   ! 61.4, 65.0 and 68.9 fail), so attempting the solve further up only spends
   ! time -- the cost grows like NMAX^4 -- on nodes that cannot be kept.
   ! Every node between X_LARGE and this ceiling must pass cross_sections_large_x.
   real(wp), parameter :: X_TM_MAX = 60.0_wp
   ! Second convergence setting.  NDGS = 4 is NOT usable: NGAUSS = NDGS*NMAX
   ! then exceeds NPNG1 = 300 well below this x ceiling, so the check would
   ! fail for a storage reason rather than a convergence one.
   real(wp), parameter :: DDELT_CHK = 3.0e-3_wp
   integer,  parameter :: NDGS_CHK  = 3
   ! Relative disagreement above which the pair is rejected.
   real(wp), parameter :: TOL_CHK   = 0.05_wp
   ! Absolute floor [Q units] on the comparison scale of the polarized
   ! channels.  Q_pol,ext and Q_bir,ext are differences of order-unity numbers,
   ! so a bare relative test on them is meaningless once they fall to the noise
   ! of the difference itself.  0.001 in Q corresponds to |C_pol|/C_ext ~ 5e-4
   ! after the size integral, i.e. ~1% of the V-band dichroism -- below the
   ! level at which this band matters at all.
   real(wp), parameter :: Q_POL_FLOOR = 1.0e-3_wp

   real(wp) :: a_eff(NA)
   real(wp), allocatable :: lambda(:), nr_cache(:), ki_cache(:)
   complex(wp), allocatable :: m_cache(:)
   real(wp) :: nr, ki

   integer, parameter :: MODE_FULL=0, MODE_TEST=1, MODE_RANGE=2, MODE_MERGE=3
   integer, parameter :: MODE_LAM=4, MODE_LAMMERGE=5
   integer  :: mode
   integer  :: jw_lo, jw_hi
   integer  :: ja_lo, ja_hi            ! a_eff window, lam modes only
   integer  :: nw                      ! wavelength-axis length of the active band
   integer  :: iarg0                   ! first argument of the mode word
   logical  :: euv_band
   ! Attempt the T-matrix above X_LARGE under the two-setting certification
   ! instead of dropping straight to geometric optics.  On for the named-
   ! wavelength modes and the EUV axis; off for the plain and `range` sweeps,
   ! which reproduce the shipped table.
   logical  :: certify_large_x
   logical  :: ja_windowed             ! this invocation covers only part of a_eff
   ! Census of the nodes offered to the two-setting certification and of those
   ! it rejected, recorded in the file header so that `lammerge` can add up the
   ! windows of one wavelength.
   integer  :: n_try_cert = 0, n_rej_cert = 0
   character(len=256) :: f_out
   character(len=256) :: f_wave_use, f_out_full_use, f_out_test_use
   character(len=256) :: lam_stem      ! output stem of a named-wavelength run
   character(len=64)  :: arg, arg2, arg3
   integer  :: ios, jw

   ! ---- CLI ---------------------------------------------------------------
   ! An optional leading `euv` switches the wavelength axis; the mode word and
   ! its arguments follow it unchanged, so the pre-existing command lines keep
   ! their exact meaning.
   euv_band = .false.
   iarg0    = 1
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg)
      if (trim(arg) == 'euv') then
         euv_band = .true.
         iarg0    = 2
      end if
   end if

   if (euv_band) then
      f_wave_use     = f_wave_euv
      f_out_full_use = f_out_full_euv
      f_out_test_use = f_out_test_euv
   else
      f_wave_use     = f_wave
      f_out_full_use = f_out_full
      f_out_test_use = f_out_test
   end if

   mode  = MODE_FULL
   ja_lo = 1
   ja_hi = NA
   ja_windowed = .false.
   certify_large_x = euv_band
   if (command_argument_count() >= iarg0) then
      call get_command_argument(iarg0, arg)
      select case (trim(arg))
      case ('test')
         mode = MODE_TEST
      case ('range')
         mode = MODE_RANGE
      case ('merge')
         mode = MODE_MERGE
      case ('lam', 'lamfile')
         mode = MODE_LAM
      case ('lammerge')
         mode = MODE_LAMMERGE
      case default
         write(*,'(a,a,a)') ' unknown mode "', trim(arg), &
              '" -- expected one of: (none), test, range, merge, lam, '// &
              'lamfile, lammerge'
         stop 1
      end select
   end if

   ! The named-wavelength modes carry their own axis, so `euv` -- which only
   ! selects between two precomputed axis files -- would have nothing to say.
   if (euv_band .and. (mode == MODE_LAM .or. mode == MODE_LAMMERGE)) then
      write(*,'(a)') ' `euv` selects a wavelength AXIS FILE and cannot be '// &
         'combined with lam/lamfile/lammerge,'
      write(*,'(a)') ' which name their wavelengths directly.'
      stop 1
   end if

   ! ---- wavelength axis ---------------------------------------------------
   ! Named-wavelength modes build it from the command line; the others count it
   ! out of the axis file so the sweep length is a property of the data.
   if (mode == MODE_LAM .or. mode == MODE_LAMMERGE) then
      certify_large_x = .true.
      call read_lam_selection()
      f_wave_use = trim(lam_stem)//'.wave'
   else
      nw = count_one_col(trim(f_wave_use))
   end if

   jw_lo = 1
   jw_hi = nw
   if (mode == MODE_RANGE) then
      if (command_argument_count() < iarg0+2) then
         write(*,'(a)') ' usage: run_q_jori.x [euv] range JW1 JW2'
         stop 1
      end if
      call get_command_argument(iarg0+1, arg2)
      call get_command_argument(iarg0+2, arg3)
      read(arg2,*,iostat=ios) jw_lo
      if (ios /= 0) then; write(*,'(a)') ' bad JW1'; stop 1; end if
      read(arg3,*,iostat=ios) jw_hi
      if (ios /= 0) then; write(*,'(a)') ' bad JW2'; stop 1; end if
      if (jw_lo < 1 .or. jw_hi > nw .or. jw_lo > jw_hi) then
         write(*,'(a,i0,a,i0,a,i0,a)') ' JW range out of bounds: [', &
            jw_lo, ', ', jw_hi, '] not in [1, ', nw, ']'
         stop 1
      end if
   end if

   ! The merge modes do not need the optics core; they only reassemble files.
   if (mode == MODE_MERGE) then
      call merge_range_files()
      stop 0
   end if
   if (mode == MODE_LAMMERGE) then
      call merge_lam_files()
      stop 0
   end if

   call read_one_col(f_aeff, NA, a_eff)
   if (mode /= MODE_LAM) then
      allocate(lambda(nw))
      call read_one_col(trim(f_wave_use), nw, lambda)
   end if
   allocate(nr_cache(nw), ki_cache(nw), m_cache(nw))
   call load_index(f_index)
   do jw = 1, nw
      call interp_m(lambda(jw), nr, ki)
      nr_cache(jw) = nr
      ki_cache(jw) = ki
      m_cache(jw)  = cmplx(nr, ki, kind=wp)
   end do

   select case (mode)
   case (MODE_TEST)
      call run_test()
   case (MODE_RANGE)
      if (euv_band) then
         write(f_out,'(a,i0,a,i0,a)') trim(f_stem_euv)//'.jw', jw_lo, '-', jw_hi, '.dat'
      else
         write(f_out,'(a,i0,a,i0,a)') trim(f_stem)//'.jw', jw_lo, '-', jw_hi, '.dat'
      end if
      call sweep_and_write(jw_lo, jw_hi, trim(f_out))
   case (MODE_LAM)
      ! The axis file is the companion the table reader needs, so it must
      ! describe the whole selection.  Only the process that owns the first
      ! a_eff column writes it, which makes concurrent windows collision-free.
      if (ja_lo == 1) call write_wave_file(trim(lam_stem)//'.wave')
      if (ja_windowed) then
         write(f_out,'(a,i0,a,i0,a)') trim(lam_stem)//'.ja', ja_lo, '-', ja_hi, '.dat'
      else
         f_out = trim(lam_stem)//'.dat'
      end if
      call sweep_and_write(1, nw, trim(f_out))
   case default
      call sweep_and_write(1, nw, trim(f_out_full_use))
   end select

contains

   subroutine sweep_and_write(jlo, jhi, fname)
      ! Full/range sweep: compute the three orientations at every (jw in
      ! [jlo,jhi], ja) node and write them in (iq, jori, jw) stream order.
      integer,          intent(in) :: jlo, jhi
      character(len=*), intent(in) :: fname
      real(wp), allocatable :: qe(:,:,:), qa(:,:,:), qs(:,:,:), qr(:,:,:)  ! (NA, NORI, jlo:jhi)
      real(wp) :: qext_ori(3), qabs_ori(3), qsca_ori(3), qre_ori(3)
      real(wp) :: x
      integer  :: ja, jjw, flag, u, n_total, n_done
      integer  :: n_try, n_rej
      logical  :: kept

      allocate(qe(NA, NORI, jlo:jhi), qa(NA, NORI, jlo:jhi), qs(NA, NORI, jlo:jhi), &
               qr(NA, NORI, jlo:jhi))
      qe = 0.0_wp;  qa = 0.0_wp;  qs = 0.0_wp;  qr = 0.0_wp

      n_total = (jhi - jlo + 1) * (ja_hi - ja_lo + 1)
      n_done  = 0
      n_try_cert = 0
      n_rej_cert = 0
      write(*,'(a)')        ' mode = sweep'
      write(*,'(a,l1)')     ' EUV band = ', euv_band
      write(*,'(a,l1)')     ' certify x > 50 against a second setting = ', certify_large_x
      write(*,'(a,a)')      ' wavelength axis = ', trim(f_wave_use)
      write(*,'(a,i0,a,i0,a)') ' jw range = [', jlo, ', ', jhi, ']'
      write(*,'(a,i0,a,i0,a)') ' ja range = [', ja_lo, ', ', ja_hi, ']'
      write(*,'(a,i0)')     ' total points = ', n_total
      write(*,'(a,a)')      ' output = ', trim(fname)

      do jjw = jlo, jhi
         n_try = 0
         n_rej = 0
         do ja = ja_lo, ja_hi
            x = 2.0_wp*PI*a_eff(ja)/lambda(jjw)
            if (certify_large_x .and. x > X_LARGE .and. x <= X_TM_MAX) then
               call cross_sections_large_x(a_eff(ja), lambda(jjw), m_cache(jjw), &
                        nr_cache(jjw), ki_cache(jjw), &
                        qext_ori, qabs_ori, qsca_ori, qre_ori, kept)
               n_try = n_try + 1
               if (.not. kept) n_rej = n_rej + 1
            else
               call oriented_cross_sections(a_eff(ja), lambda(jjw), m_cache(jjw), &
                        EPS_BA, NP_OBL, DDELT, NDGS, qext_ori, qabs_ori, qsca_ori, flag, &
                        qre_ori=qre_ori)
            end if
            qe(ja, :, jjw) = qext_ori
            qa(ja, :, jjw) = qabs_ori
            qs(ja, :, jjw) = qsca_ori
            qr(ja, :, jjw) = qre_ori
            n_done = n_done + 1
         end do
         n_try_cert = n_try_cert + n_try
         n_rej_cert = n_rej_cert + n_rej
         if (certify_large_x) &
            write(*,'(a,i0,a,es12.5,a,i0,a,i0)') ' jw ', jjw, '  lam=', lambda(jjw), &
               ' um   x>50 attempted = ', n_try, '   rejected = ', n_rej
         if (mod(jjw - jlo + 1, max((jhi-jlo+1)/20, 1)) == 0) then
            write(*,'(a,i0,a,i0,a,f6.1,a)') ' progress: ', n_done, '/', &
               n_total, '  (', 100.0_wp*real(n_done,wp)/real(n_total,wp), '%)'
         end if
      end do

      open(newunit=u, file=trim(fname), status='replace', action='write')
      call write_header(u, jlo, jhi)
      call write_block(u, qe, jlo, jhi)   ! iq = 1  Q_ext
      call write_block(u, qa, jlo, jhi)   ! iq = 2  Q_abs
      call write_block(u, qs, jlo, jhi)   ! iq = 3  Q_sca
      call write_block(u, qr, jlo, jhi)   ! iq = 4  Q_re  (birefringence twin)
      close(u)
      write(*,'(a,a)') ' wrote ', trim(fname)
      if (certify_large_x) &
         write(*,'(a,i0,a,i0)') ' LARGE_X_CENSUS attempted = ', n_try_cert, &
            '   rejected = ', n_rej_cert

      deallocate(qe, qa, qs, qr)
   end subroutine sweep_and_write


   subroutine cross_sections_large_x(a_um, lam_um, m_ref, n_r, k_i, &
                                     qext_ori, qabs_ori, qsca_ori, qre_ori, kept)
      ! Orientation-resolved optics at a size parameter beyond x = 50, where
      ! the T-matrix is no longer trusted on its own status flag.
      !
      ! WHY.  Mishchenko's solver reports IERR = 0 once its own convergence
      ! test is satisfied at the multipole order it reached, but the storage
      ! limit NPN1 = 100 bounds that order, so for x >~ 55 the test can be met
      ! by an expansion that is truncated rather than converged.  A measured
      ! instance: (lambda, a_eff) = (0.0912, 0.9441) um, x = 65.0, returns
      ! IERR = 0 with Q_pol,ext = -1.16e-2 while its x = 58 and x = 61
      ! neighbours give -6.5e-3, and a different tolerance turns the same node
      ! into IERR = 3.  A 50%-wrong number therefore arrives labelled
      ! converged.
      !
      ! WHAT IS DONE.  The node is solved twice under independent convergence
      ! settings -- the production (DDELT, NDGS) and a coarser-quadrature,
      ! looser-tolerance pair -- and the T-matrix answer is kept only if the
      ! two agree to TOL_CHK on every channel: the three Q_ext, the three
      ! Q_abs, the three Q_sca, and the polarized combinations
      ! Q_pol,ext = 0.5*(Q_ext(3)-Q_ext(2)) and Q_bir,ext = 0.5*(Q_re(3)-Q_re(2)).
      ! Two truncated expansions of different order do not agree; two converged
      ! ones do.  On disagreement (or on either solve failing) the node falls
      ! to geometric_optics_limit, the correct asymptote for x >> 1, which
      ! carries Q_pol,ext = Q_bir,ext = 0 exactly -- a documented deficit rather
      ! than an unconverged number.
      !
      ! NO EXTRAPOLATION is performed above the certified region: the missing
      ! dichroic extinction is left at the geometric-optics zero.
      real(wp),    intent(in)  :: a_um, lam_um, n_r, k_i
      complex(wp), intent(in)  :: m_ref
      real(wp),    intent(out) :: qext_ori(3), qabs_ori(3), qsca_ori(3), qre_ori(3)
      logical,     intent(out) :: kept

      real(wp) :: qe1(3), qa1(3), qs1(3), qr1(3)
      real(wp) :: qe2(3), qa2(3), qs2(3), qr2(3)
      real(wp) :: qext_go, qsca_go, walb_go, asymm_go
      integer  :: ierr1, ierr2, j

      kept = .false.
      call tmatrix_oriented_cross(a_um, lam_um, m_ref, EPS_BA, NP_OBL, DDELT, NDGS, &
                                  qe1, qs1, qa1, ierr1, qre_ori=qr1)
      if (ierr1 == 0) then
         call tmatrix_oriented_cross(a_um, lam_um, m_ref, EPS_BA, NP_OBL, &
                                     DDELT_CHK, NDGS_CHK, qe2, qs2, qa2, ierr2, &
                                     qre_ori=qr2)
         if (ierr2 == 0) then
            kept = .true.
            ! An identically zero extinction means nothing was reconstructed.
            if (maxval(abs(qe1)) <= 0.0_wp) kept = .false.
            do j = 1, 3
               if (.not. agree(qe1(j), qe2(j), TOL_CHK, 0.0_wp)) kept = .false.
               if (.not. agree(qa1(j), qa2(j), TOL_CHK, 0.0_wp)) kept = .false.
               if (.not. agree(qs1(j), qs2(j), TOL_CHK, 0.0_wp)) kept = .false.
            end do
            if (.not. agree(0.5_wp*(qe1(3)-qe1(2)), 0.5_wp*(qe2(3)-qe2(2)), &
                            TOL_CHK, Q_POL_FLOOR)) kept = .false.
            if (.not. agree(0.5_wp*(qr1(3)-qr1(2)), 0.5_wp*(qr2(3)-qr2(2)), &
                            TOL_CHK, Q_POL_FLOOR)) kept = .false.
            do j = 1, 3
               if (.not. ieee_is_finite(qe1(j))) kept = .false.
               if (.not. ieee_is_finite(qa1(j))) kept = .false.
               if (.not. ieee_is_finite(qs1(j))) kept = .false.
               if (.not. ieee_is_finite(qr1(j))) kept = .false.
            end do
         end if
      end if

      if (kept) then
         qext_ori = qe1;  qabs_ori = qa1;  qsca_ori = qs1;  qre_ori = qr1
      else
         call geometric_optics_limit(a_um, lam_um, n_r, k_i, EPS_BA, &
                 qext_go, qsca_go, walb_go, asymm_go, &
                 qext_ori=qext_ori, qabs_ori=qabs_ori, qsca_ori=qsca_ori, &
                 qre_ori=qre_ori)
      end if
   end subroutine cross_sections_large_x


   logical function agree(v1, v2, tol, floor_abs)
      ! |v1-v2| measured against max(|v1|, |v2|, floor_abs).  floor_abs = 0
      ! makes it a plain relative test; a positive floor stops a difference
      ! that is already negligible in absolute terms from failing on a
      ! meaningless relative scale.
      real(wp), intent(in) :: v1, v2, tol, floor_abs
      real(wp) :: s
      s = max(abs(v1), abs(v2), floor_abs, tiny(1.0_wp))
      agree = (abs(v1 - v2) <= tol * s)
   end function agree


   subroutine write_block(u, q, jlo, jhi)
      ! One quantity: for jori = 1..3, for jw = jlo..jhi, one record of the
      ! a_eff columns this invocation owns (all NA of them unless ja= narrowed
      ! it).  Matches the inner two loops of q_table_jori.f90:373-396.
      integer,  intent(in) :: u, jlo, jhi
      real(wp), intent(in) :: q(NA, NORI, jlo:jhi)
      integer :: jori, jjw
      do jori = 1, NORI
         do jjw = jlo, jhi
            write(u,'(*(es13.5))') q(ja_lo:ja_hi, jori, jjw)
         end do
      end do
   end subroutine write_block


   subroutine write_header(u, jlo, jhi)
      ! Exactly 12 header lines, so load_q_table_jori (NHEAD = 12) skips them.
      ! Line 11 names the wavelengths of a named-wavelength product, and line 12
      ! carries the machine-readable jw and a_eff windows plus the large-x
      ! census, which the merge modes parse.
      integer, intent(in) :: u, jlo, jhi
      write(u,'(a)') '# DH21 astrodust orientation-resolved Q, P = 0.20, fFe = 0.00, b/a = 1.4'
      write(u,'(a,i0,a,i0,a)') '# Q = C/(pi a_eff^2); a_eff from DH21_aeff (', NA, &
         '), lambda from '//trim(basename(f_wave_use))//' (', nw, ').'
      write(u,'(a)') '# jori convention (a = spheroid symmetry axis):'
      write(u,'(a)') '#   jori=1: k || a'
      write(u,'(a)') '#   jori=2: k perp a, E || a'
      write(u,'(a)') '#   jori=3: k perp a, E perp a'
      write(u,'(a,i0,a)') '# Stream order after this header (free format, one record = ', &
         ja_hi - ja_lo + 1, ' a_eff values):'
      write(u,'(a,i0,a)') '#   do iq=1,4 ; do jori=1,3 ; do jw=1,', nw, ' ; write one record ; end'
      write(u,'(a)') '#   iq = 1 Q_ext, 2 Q_abs, 3 Q_sca, 4 Q_re (birefringence twin)'
      if (certify_large_x) then
         write(u,'(a,f5.1,a)') '# Regime: x<0.1 Rayleigh, 0.1<=x<=50 T-matrix, '// &
            '50<x<=', X_TM_MAX, ' T-matrix if two-setting certified, else geometric optics.'
      else
         write(u,'(a)') '# Regime: x<0.1 Rayleigh, 0.1<=x<=50 T-matrix, x>50 geometric optics.'
      end if
      if (mode == MODE_LAM .or. mode == MODE_LAMMERGE) then
         write(u,'(a,*(1x,es20.12))') '# LAMBDA_UM', lambda(1:nw)
      else
         write(u,'(a)') '# Generated by tmatrix/driver/run_q_jori.f90.'
      end if
      write(u,'(a,i0,1x,i0,a,i0,1x,i0,a,i0,1x,i0)') '# JW_WINDOW ', jlo, jhi, &
         '  JA_WINDOW ', ja_lo, ja_hi, '  LARGE_X_CENSUS ', n_try_cert, n_rej_cert
   end subroutine write_header


   function basename(path) result(b)
      ! Trailing path component, for the header line that names the axis file.
      character(len=*), intent(in) :: path
      character(len=len(path))     :: b
      integer :: k
      k = index(trim(path), '/', back=.true.)
      b = trim(path(k+1:))
   end function basename


   subroutine run_test()
      ! Strided ~7x7 sample.  Times each node, classifies it as cheap
      ! (Rayleigh / geometric optics) or expensive (T-matrix attempted), and
      ! extrapolates the two averages onto the full grid for a wall-time
      ! estimate.  Writes a diagnostic columnar file (not the stream format).
      integer  :: ja_step, jw_step
      integer  :: ja, jjw, flag, u
      real(wp) :: qext_ori(3), qabs_ori(3), qsca_ori(3), qre_ori(3), x
      integer(int64) :: c0, c1, crate
      real(wp) :: dt, x_tm_hi
      integer  :: n_cheap_s, n_exp_s
      real(wp) :: t_cheap_s, t_exp_s
      integer  :: n_cheap_f, n_exp_f
      real(wp) :: avg_cheap, avg_exp, est
      logical  :: expensive, kept

      ja_step = 28
      jw_step = max(nw/6, 1)
      x_tm_hi = X_LARGE
      if (certify_large_x) x_tm_hi = X_TM_MAX

      call system_clock(count_rate=crate)

      ! Full-grid regime census (arithmetic only, no solve).
      n_cheap_f = 0;  n_exp_f = 0
      do jjw = 1, nw
         do ja = 1, NA
            x = 2.0_wp*PI*a_eff(ja)/lambda(jjw)
            if (x < X_SMALL .or. x > x_tm_hi) then
               n_cheap_f = n_cheap_f + 1
            else
               n_exp_f = n_exp_f + 1
            end if
         end do
      end do

      open(newunit=u, file=trim(f_out_test_use), status='replace', action='write')
      write(u,'(a)') '# run_q_jori test sample (strided).  Columns:'
      write(u,'(a)') '#  lambda[um]  a_eff[um]  x  flag  Qext(1:3)  Qabs(1:3)  Qsca(1:3)'

      n_cheap_s = 0;  n_exp_s = 0
      t_cheap_s = 0.0_wp;  t_exp_s = 0.0_wp
      write(*,'(a)')    ' mode = test (strided sample)'
      write(*,'(a,l1)') ' EUV band = ', euv_band
      do jjw = 1, nw, jw_step
         do ja = 1, NA, ja_step
            x = 2.0_wp*PI*a_eff(ja)/lambda(jjw)
            expensive = (x >= X_SMALL .and. x <= x_tm_hi)
            call system_clock(c0)
            if (certify_large_x .and. x > X_LARGE .and. x <= X_TM_MAX) then
               call cross_sections_large_x(a_eff(ja), lambda(jjw), m_cache(jjw), &
                        nr_cache(jjw), ki_cache(jjw), &
                        qext_ori, qabs_ori, qsca_ori, qre_ori, kept)
               flag = 0;  if (.not. kept) flag = 20
            else
               call oriented_cross_sections(a_eff(ja), lambda(jjw), m_cache(jjw), &
                        EPS_BA, NP_OBL, DDELT, NDGS, qext_ori, qabs_ori, qsca_ori, flag)
            end if
            call system_clock(c1)
            dt = real(c1 - c0, kind=wp) / real(crate, kind=wp)
            if (expensive) then
               n_exp_s = n_exp_s + 1;  t_exp_s = t_exp_s + dt
            else
               n_cheap_s = n_cheap_s + 1;  t_cheap_s = t_cheap_s + dt
            end if
            write(u,'(2es13.5,es12.4,i5,9es13.5)') lambda(jjw), a_eff(ja), x, flag, &
               qext_ori, qabs_ori, qsca_ori
            write(*,'(a,es11.4,a,es11.4,a,es10.3,a,i4,a,f8.2,a)') &
               '  lam=', lambda(jjw), ' a=', a_eff(ja), ' x=', x, ' flag=', flag, &
               '  (', dt, ' s)'
         end do
      end do
      close(u)

      avg_cheap = 0.0_wp;  if (n_cheap_s > 0) avg_cheap = t_cheap_s/real(n_cheap_s,wp)
      avg_exp   = 0.0_wp;  if (n_exp_s   > 0) avg_exp   = t_exp_s  /real(n_exp_s,wp)
      est = real(n_cheap_f,wp)*avg_cheap + real(n_exp_f,wp)*avg_exp

      write(*,'(a)')            ' --- full-sweep wall-time estimate ---'
      write(*,'(a,a)')          ' wrote ', trim(f_out_test_use)
      write(*,'(a,i0,a,es10.3,a)') ' cheap nodes sampled : ', n_cheap_s, &
         '   avg ', avg_cheap, ' s'
      write(*,'(a,i0,a,es10.3,a)') ' T-matrix nodes samp.: ', n_exp_s, &
         '   avg ', avg_exp, ' s'
      write(*,'(a,i0,a,i0)')    ' full grid: cheap = ', n_cheap_f, '   T-matrix = ', n_exp_f
      write(*,'(a,f10.1,a,f8.2,a,f7.3,a)') ' estimated full sweep: ', est, &
         ' s  = ', est/60.0_wp, ' min  = ', est/3600.0_wp, ' h'
   end subroutine run_test


   subroutine merge_range_files()
      ! Reads the range files named on the command line (after the mode word),
      ! places each jw record into the full (iq, jori, jw) arrays, verifies the
      ! union covers jw = 1..nw with no gaps or overlaps, and writes the full
      ! file in q_table_jori.f90 stream order.
      real(wp), allocatable :: qe(:,:,:), qa(:,:,:), qs(:,:,:), qr(:,:,:)  ! (NA, NORI, nw)
      logical, allocatable  :: covered(:)
      integer  :: nfiles, k, u, ios, i, iq, jori, jjw, jlo, jhi
      character(len=256) :: path
      character(len=512) :: line
      real(wp) :: row(NA)

      nfiles = command_argument_count() - iarg0
      if (nfiles < 1) then
         write(*,'(a)') ' usage: run_q_jori.x [euv] merge FILE [FILE ...]'
         stop 1
      end if

      allocate(qe(NA, NORI, nw), qa(NA, NORI, nw), qs(NA, NORI, nw), qr(NA, NORI, nw))
      allocate(covered(nw))
      qe = 0.0_wp;  qa = 0.0_wp;  qs = 0.0_wp;  qr = 0.0_wp
      covered = .false.

      do k = 1, nfiles
         call get_command_argument(iarg0+k, path)
         open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
         if (ios /= 0) then
            write(error_unit,'(a,a)') ' merge: cannot open ', trim(path)
            stop 1
         end if
         ! Read 12 header lines; recover the jw window from the JW_WINDOW line.
         jlo = 0;  jhi = -1
         do i = 1, 12
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) then
               write(error_unit,'(a,a)') ' merge: short header in ', trim(path)
               stop 1
            end if
            if (index(line, 'JW_WINDOW') > 0) &
               read(line(index(line,'JW_WINDOW')+9:), *) jlo, jhi
         end do
         if (jlo < 1 .or. jhi > nw .or. jlo > jhi) then
            write(error_unit,'(a,a)') ' merge: bad JW_WINDOW in ', trim(path)
            stop 1
         end if
         ! Body: (iq, jori, jw in [jlo,jhi]) records of NA values.
         do iq = 1, 4
            do jori = 1, NORI
               do jjw = jlo, jhi
                  read(u,*,iostat=ios) row
                  if (ios /= 0) then
                     write(error_unit,'(a,a)') ' merge: short body in ', trim(path)
                     stop 1
                  end if
                  select case (iq)
                  case (1);  qe(:, jori, jjw) = row
                  case (2);  qa(:, jori, jjw) = row
                  case (3);  qs(:, jori, jjw) = row
                  case (4);  qr(:, jori, jjw) = row
                  end select
               end do
            end do
         end do
         close(u)
         do jjw = jlo, jhi
            if (covered(jjw)) then
               write(error_unit,'(a,i0)') ' merge: jw covered twice at jw=', jjw
               stop 1
            end if
            covered(jjw) = .true.
         end do
         write(*,'(a,a,a,i0,a,i0,a)') ' read ', trim(path), '  (jw ', jlo, '..', jhi, ')'
      end do

      do jjw = 1, nw
         if (.not. covered(jjw)) then
            write(error_unit,'(a,i0)') ' merge: jw not covered at jw=', jjw
            stop 1
         end if
      end do

      open(newunit=u, file=trim(f_out_full_use), status='replace', action='write')
      call write_header(u, 1, nw)
      call write_block(u, qe, 1, nw)
      call write_block(u, qa, 1, nw)
      call write_block(u, qs, 1, nw)
      call write_block(u, qr, 1, nw)
      close(u)
      write(*,'(a,a)') ' wrote ', trim(f_out_full_use)
      deallocate(qe, qa, qs, qr, covered)
   end subroutine merge_range_files


   subroutine read_lam_selection()
      ! Builds the wavelength list of a named-wavelength run, the a_eff window
      ! this invocation owns, and the stem its products are written under.
      !
      !   lam       the wavelengths [um] are the remaining arguments
      !   lamfile   they are the numbers in the one file named
      !   lammerge  the first remaining argument is the output STEM, and the
      !             wavelengths come from its companion STEM.wave, so a merge
      !             never has to be told the axis a second time
      !
      ! ja=JA1:JA2 and tag=NAME may appear anywhere among them.  The list is
      ! sorted ascending and must be free of duplicates, because the table
      ! reader requires a strictly increasing axis.
      character(len=256) :: mode_word, a, listfile, tag
      real(wp), allocatable :: vals(:)
      real(wp) :: v
      integer  :: k, nargs, nval, ios2
      logical  :: from_file

      call get_command_argument(iarg0, mode_word)
      from_file = (trim(mode_word) == 'lamfile')
      nargs     = command_argument_count()
      tag       = ''
      listfile  = ''
      lam_stem  = ''
      allocate(vals(max(nargs, 1)))
      nval = 0

      do k = iarg0+1, nargs
         call get_command_argument(k, a)
         if (a(1:3) == 'ja=') then
            call parse_ja_window(trim(a(4:)))
         else if (a(1:4) == 'tag=') then
            tag = trim(a(5:))
         else if (mode == MODE_LAMMERGE) then
            ! The first is the stem; the rest are the files merge_lam_files reads.
            if (len_trim(lam_stem) == 0) lam_stem = trim(a)
         else if (from_file) then
            if (len_trim(listfile) /= 0) then
               write(*,'(a,a)') ' lamfile takes ONE path; unexpected argument ', trim(a)
               stop 1
            end if
            listfile = trim(a)
         else
            read(a,*,iostat=ios2) v
            if (ios2 /= 0 .or. v <= 0.0_wp) then
               write(*,'(a,a,a)') ' "', trim(a), '" is not a positive wavelength [um]'
               stop 1
            end if
            nval = nval + 1
            vals(nval) = v
         end if
      end do

      if (mode == MODE_LAMMERGE) then
         if (len_trim(lam_stem) == 0) then
            write(*,'(a)') ' usage: run_q_jori.x lammerge STEM FILE [FILE ...]'
            stop 1
         end if
         nw = count_one_col(trim(lam_stem)//'.wave')
         allocate(lambda(nw))
         call read_one_col(trim(lam_stem)//'.wave', nw, lambda)
         deallocate(vals)
         return
      end if

      if (from_file) then
         if (len_trim(listfile) == 0) then
            write(*,'(a)') ' usage: run_q_jori.x lamfile PATH [ja=JA1:JA2] [tag=NAME]'
            stop 1
         end if
         deallocate(vals)
         call read_lam_list(trim(listfile), vals, nval)
      end if

      if (nval < 1) then
         write(*,'(a)') ' no wavelengths given.'
         write(*,'(a)') ' usage: run_q_jori.x lam L1 [L2 ...] [ja=JA1:JA2] [tag=NAME]'
         stop 1
      end if

      call sort_ascending(vals, nval)
      do k = 2, nval
         if (vals(k) <= vals(k-1)) then
            write(*,'(a,es13.5,a)') ' duplicate wavelength ', vals(k), &
               ' um: the axis must be strictly increasing.'
            stop 1
         end if
      end do

      nw = nval
      allocate(lambda(nw))
      lambda = vals(1:nw)
      deallocate(vals)

      if (len_trim(tag) == 0) then
         if (nw == 1) then
            tag = 'lam'//lam_token(lambda(1))
         else
            write(a,'(i0)') nw
            tag = 'lam'//lam_token(lambda(1))//'_'//lam_token(lambda(nw))// &
                  '_n'//trim(a)
         end if
      end if
      lam_stem = f_stem//'.'//trim(tag)

      write(*,'(a,i0,a)') ' named-wavelength run: ', nw, ' wavelength(s) [um]'
      do k = 1, nw
         write(*,'(a,i0,a,es14.6)') '   jw ', k, ' : ', lambda(k)
      end do
      write(*,'(a,a)') ' output stem = ', trim(lam_stem)
   end subroutine read_lam_selection


   subroutine parse_ja_window(s)
      ! ja=JA1:JA2 -- the a_eff columns this process computes.  Every
      ! (lambda, a_eff) node is solved on its own, so splitting the 169 radii
      ! over processes changes the wall time and nothing else.
      character(len=*), intent(in) :: s
      integer :: p, ios2
      p = index(s, ':')
      if (p < 2 .or. p >= len_trim(s)) then
         write(*,'(a,a)') ' expected ja=JA1:JA2, got ja=', trim(s)
         stop 1
      end if
      read(s(1:p-1),*,iostat=ios2) ja_lo
      if (ios2 /= 0) then; write(*,'(a)') ' bad JA1'; stop 1; end if
      read(s(p+1:),*,iostat=ios2) ja_hi
      if (ios2 /= 0) then; write(*,'(a)') ' bad JA2'; stop 1; end if
      if (ja_lo < 1 .or. ja_hi > NA .or. ja_lo > ja_hi) then
         write(*,'(a,i0,a,i0,a,i0,a)') ' JA window out of bounds: [', ja_lo, &
            ', ', ja_hi, '] not in [1, ', NA, ']'
         stop 1
      end if
      ja_windowed = .true.
   end subroutine parse_ja_window


   subroutine read_lam_list(path, v, n)
      ! Wavelengths [um] from a text file: '#' starts a comment, blank lines
      ! are ignored, and any number of values may share a line (up to 4096
      ! characters of it).  Counted on a first pass, read on a second.
      character(len=*),      intent(in)  :: path
      real(wp), allocatable, intent(out) :: v(:)
      integer,               intent(out) :: n
      integer :: u, ios2, ipass, k, ic, nt
      character(len=4096) :: buf

      do ipass = 1, 2
         open(newunit=u, file=path, status='old', action='read', iostat=ios2)
         if (ios2 /= 0) then
            write(*,'(a,a)') ' ERROR: cannot open ', trim(path)
            stop 1
         end if
         n = 0
         do
            read(u,'(a)',iostat=ios2) buf
            if (ios2 /= 0) exit
            ic = index(buf, '#')
            if (ic > 0) buf(ic:) = ' '
            nt = count_tokens(buf)
            if (nt == 0) cycle
            if (ipass == 2) then
               read(buf,*,iostat=ios2) (v(k), k = n+1, n+nt)
               if (ios2 /= 0) then
                  write(*,'(a,a)') ' ERROR: non-numeric wavelength in ', trim(path)
                  stop 1
               end if
            end if
            n = n + nt
         end do
         close(u)
         if (ipass == 1) allocate(v(max(n,1)))
      end do
   end subroutine read_lam_list


   integer function count_tokens(s) result(n)
      ! Whitespace-separated tokens in one string.
      character(len=*), intent(in) :: s
      integer :: i
      logical :: in_tok
      n = 0
      in_tok = .false.
      do i = 1, len_trim(s)
         if (s(i:i) == ' ' .or. s(i:i) == char(9)) then
            in_tok = .false.
         else
            if (.not. in_tok) n = n + 1
            in_tok = .true.
         end if
      end do
   end function count_tokens


   subroutine sort_ascending(v, n)
      ! Insertion sort; n is a handful of wavelengths at most.
      real(wp), intent(inout) :: v(:)
      integer,  intent(in)    :: n
      integer  :: i, k
      real(wp) :: t
      do i = 2, n
         t = v(i)
         k = i - 1
         do while (k >= 1)
            if (v(k) <= t) exit
            v(k+1) = v(k)
            k = k - 1
         end do
         v(k+1) = t
      end do
   end subroutine sort_ascending


   function lam_token(x) result(s)
      ! Wavelength as it appears in a file name: 9 characters, e.g. 6.020E-02.
      real(wp), intent(in) :: x
      character(len=9) :: s
      write(s,'(es9.3)') x
   end function lam_token


   subroutine write_wave_file(path)
      ! The companion wavelength axis of a named-wavelength product, in the
      ! format read_grid / count_grid_values expect (two title lines, then the
      ! values).  Written 8 per line: a single very long record is legal
      ! Fortran but is what a fixed-buffer token count truncates.
      character(len=*), intent(in) :: path
      integer :: u, k
      open(newunit=u, file=path, status='replace', action='write')
      write(u,'(a)') 'DH21Ad...'
      write(u,'(a,i0,a)') 'List of wavelengths(um) wave(jw), jw=0-', nw-1, &
         '  [named nodes, written by run_q_jori.f90]'
      write(u,'(8es20.12)') (lambda(k), k = 1, nw)
      close(u)
      write(*,'(a,a)') ' wrote ', trim(path)
   end subroutine write_wave_file


   subroutine merge_lam_files()
      ! Assembles the a_eff windows of one named-wavelength selection.  Each
      ! input file holds every jw record but only its own a_eff columns, so the
      ! windows are placed side by side and their union must cover 1..NA once.
      real(wp), allocatable :: qe(:,:,:), qa(:,:,:), qs(:,:,:), qr(:,:,:)
      logical, allocatable  :: covered(:)
      integer  :: nfiles, k, u, ios2, i, iq, jori, jjw, jlo, jhi, jalo, jahi, nc
      integer  :: n_try_f, n_rej_f
      character(len=256) :: path
      character(len=512) :: line
      real(wp) :: row(NA)

      nfiles = 0
      do k = iarg0+1, command_argument_count()
         call get_command_argument(k, path)
         if (path(1:3) == 'ja=' .or. path(1:4) == 'tag=') cycle
         nfiles = nfiles + 1
      end do
      nfiles = nfiles - 1                       ! the first non-option is the stem
      if (nfiles < 1) then
         write(*,'(a)') ' usage: run_q_jori.x lammerge STEM FILE [FILE ...]'
         stop 1
      end if

      allocate(qe(NA, NORI, nw), qa(NA, NORI, nw), qs(NA, NORI, nw), qr(NA, NORI, nw))
      allocate(covered(NA))
      qe = 0.0_wp;  qa = 0.0_wp;  qs = 0.0_wp;  qr = 0.0_wp
      covered = .false.
      n_try_cert = 0
      n_rej_cert = 0

      nc = 0
      do k = iarg0+1, command_argument_count()
         call get_command_argument(k, path)
         if (path(1:3) == 'ja=' .or. path(1:4) == 'tag=') cycle
         nc = nc + 1
         if (nc == 1) cycle                     ! the stem
         open(newunit=u, file=trim(path), status='old', action='read', iostat=ios2)
         if (ios2 /= 0) then
            write(error_unit,'(a,a)') ' lammerge: cannot open ', trim(path)
            stop 1
         end if
         jlo = 0;  jhi = -1;  jalo = 0;  jahi = -1
         n_try_f = 0;  n_rej_f = 0
         do i = 1, 12
            read(u,'(a)',iostat=ios2) line
            if (ios2 /= 0) then
               write(error_unit,'(a,a)') ' lammerge: short header in ', trim(path)
               stop 1
            end if
            if (index(line, 'JW_WINDOW') > 0) &
               read(line(index(line,'JW_WINDOW')+9:), *) jlo, jhi
            if (index(line, 'JA_WINDOW') > 0) &
               read(line(index(line,'JA_WINDOW')+9:), *) jalo, jahi
            if (index(line, 'LARGE_X_CENSUS') > 0) &
               read(line(index(line,'LARGE_X_CENSUS')+14:), *) n_try_f, n_rej_f
         end do
         if (jlo /= 1 .or. jhi /= nw) then
            write(error_unit,'(a,a)') ' lammerge: wavelength count differs from '// &
               trim(lam_stem)//'.wave in ', trim(path)
            stop 1
         end if
         if (jalo < 1 .or. jahi > NA .or. jalo > jahi) then
            write(error_unit,'(a,a)') ' lammerge: bad JA_WINDOW in ', trim(path)
            stop 1
         end if
         do iq = 1, 4
            do jori = 1, NORI
               do jjw = 1, nw
                  read(u,*,iostat=ios2) row(jalo:jahi)
                  if (ios2 /= 0) then
                     write(error_unit,'(a,a)') ' lammerge: short body in ', trim(path)
                     stop 1
                  end if
                  select case (iq)
                  case (1);  qe(jalo:jahi, jori, jjw) = row(jalo:jahi)
                  case (2);  qa(jalo:jahi, jori, jjw) = row(jalo:jahi)
                  case (3);  qs(jalo:jahi, jori, jjw) = row(jalo:jahi)
                  case (4);  qr(jalo:jahi, jori, jjw) = row(jalo:jahi)
                  end select
               end do
            end do
         end do
         close(u)
         do i = jalo, jahi
            if (covered(i)) then
               write(error_unit,'(a,i0)') ' lammerge: a_eff covered twice at ja=', i
               stop 1
            end if
            covered(i) = .true.
         end do
         n_try_cert = n_try_cert + n_try_f
         n_rej_cert = n_rej_cert + n_rej_f
         write(*,'(a,a,a,i0,a,i0,a)') ' read ', trim(path), '  (ja ', jalo, '..', jahi, ')'
      end do

      do i = 1, NA
         if (.not. covered(i)) then
            write(error_unit,'(a,i0)') ' lammerge: a_eff not covered at ja=', i
            stop 1
         end if
      end do

      ja_lo = 1;  ja_hi = NA
      open(newunit=u, file=trim(lam_stem)//'.dat', status='replace', action='write')
      call write_header(u, 1, nw)
      call write_block(u, qe, 1, nw)
      call write_block(u, qa, 1, nw)
      call write_block(u, qs, 1, nw)
      call write_block(u, qr, 1, nw)
      close(u)
      write(*,'(a,a)') ' wrote ', trim(lam_stem)//'.dat'
      write(*,'(a,i0,a,i0)') ' LARGE_X_CENSUS attempted = ', n_try_cert, &
         '   rejected = ', n_rej_cert
      deallocate(qe, qa, qs, qr, covered)
   end subroutine merge_lam_files


   subroutine read_one_col(filename, n, x)
      ! DH21_aeff / DH21_wave / DH21_wave_euv: 2 header lines, then all n
      ! values in free format.
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: n
      real(wp),         intent(out) :: x(n)
      integer :: u, ios
      character(len=512) :: header
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a,a)') ' ERROR: cannot open ', trim(filename)
         stop 1
      end if
      read(u,'(a)') header
      read(u,'(a)') header
      read(u,*) x(1:n)
      close(u)
   end subroutine read_one_col


   integer function count_one_col(filename) result(n)
      ! Number of whitespace-separated values after the 2 header lines, so the
      ! wavelength axis sets the sweep length instead of a compiled-in
      ! constant.  Used for both axes; DH21_wave returns 1129.
      !
      ! The record is consumed in fixed-size chunks with a non-advancing read
      ! rather than in one plain read, because a record longer than the buffer
      ! would otherwise be silently truncated -- and DH21_wave puts all 1129
      ! values on ONE 11289-character line, which a 4096-character buffer cut
      ! down to 410.  A token straddling a chunk boundary is not double
      ! counted: in_tok carries across chunks and is cleared only at end of
      ! record.
      character(len=*), intent(in) :: filename
      integer :: u, ios, i, L
      character(len=4096) :: line
      logical :: in_tok
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a,a)') ' ERROR: cannot open ', trim(filename)
         stop 1
      end if
      read(u,'(a)',iostat=ios) line
      read(u,'(a)',iostat=ios) line
      n = 0
      in_tok = .false.
      do
         read(u,'(a)',advance='no',size=L,iostat=ios) line
         do i = 1, L
            if (line(i:i) == ' ' .or. line(i:i) == char(9)) then
               in_tok = .false.
            else
               if (.not. in_tok) n = n + 1
               in_tok = .true.
            end if
         end do
         if (is_iostat_eor(ios)) then
            in_tok = .false.            ! end of record ends any open token
         else if (ios /= 0) then
            exit                        ! end of file
         end if
      end do
      close(u)
      if (n < 2) then
         write(*,'(a,a)') ' ERROR: no wavelength values in ', trim(filename)
         stop 1
      end if
   end function count_one_col

end program run_q_jori
