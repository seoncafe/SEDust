program spheroid_asymmetry_table
   ! Random-orientation scattering asymmetry g = <cos theta> of the two
   ! spheroid populations of the G18 Model D dust model (Guillet et al. 2018,
   ! A&A 610, A16, their model D), written in the DustEM G_<gtype>.DAT layout.
   !
   ! WHY THIS EXISTS.  The DustEM distribution ships Q_ tables for both
   ! populations but no G_ table: Guillet et al. never computed an asymmetry
   ! parameter (the paper contains no such quantity), and DustEM reads a G_
   ! file only under its `pdr` run keyword, which this model does not set.  A
   ! radiative-transfer host that is handed g = 0 for a model whose albedo
   ! reaches 0.38 scatters isotropically, which is wrong for grains that reach
   ! x = 2 pi a / lambda ~ 100 in the ultraviolet.  This program computes g
   ! with the same physics the published Q tables are built from, so that the
   ! model's scattering asymmetry is a computed quantity rather than an
   ! absence.
   !
   ! WHAT IS REPRODUCED FROM THE PUBLISHED MODEL, item by item:
   !
   !   population              material                          shape   sizes
   !   amCBE_0.3333x           BE amorphous carbon, Zubko et al. eps=1/3  100,
   !                           (1996), burning of benzene in air          2e-4
   !                                                                      to 2 um
   !   aSil2001BE6pctG_0.4x    WD01 astrosilicate matrix (the     eps=0.4 70,
   !                           pre-2001 "smoothed UV" Draine &            3e-3
   !                           Lee 1984 / Laor & Draine 1993              to 2 um
   !                           index, beta = 2) with 6% by
   !                           volume of BE a-C inclusions,
   !                           Maxwell Garnett (their Eq. 28)
   !
   !   * eps = b/a < 1 is a PROLATE spheroid in Mishchenko's convention, which
   !     is what spheroid_optics and the T-matrix library use.  Guillet et al.
   !     Table 3 row D gives 1/3 for the carbon and 1/2.5 for the composite
   !     silicate ("an elongation of 2.5 is enough", their Sect. 5.6).
   !   * Q and g are normalized to and defined against the VOLUME-EQUIVALENT
   !     sphere radius a_V (their Sect. 3.6).  The T-matrix core is called with
   !     RAT = 1, so its a_um argument is a_V and its efficiencies carry the
   !     same normalization.
   !   * The size grid is taken from the population's own Q_<gtype>.DAT so the
   !     two tables index alike, and the wavelength grid is oprop/LAMBDA.DAT.
   !   * Guillet et al. formed their orientation average as a sum over 768
   !     HEALPix directions of fixed-orientation cross sections (their
   !     Eq. 19).  Mishchenko's random-orientation T-matrix, which this
   !     program calls, is the exact analytic form of that average.
   !
   ! REGIME POLICY.  For each (a, lambda) with x = 2 pi a_V / lambda:
   !
   !   regime 1, x < 0.1
   !       g = 0 exactly.  In the Rayleigh dipole limit the scattered
   !       intensity is symmetric about the scattering plane's equator and
   !       <cos theta> vanishes; the leading correction is O(x^2).  The
   !       refractive index is not even needed, and is not read.
   !
   !   regime 2, 0.1 <= x <= 50 and the T-matrix converges
   !       g = the random-orientation T-matrix asymmetry.  This is the answer
   !       wherever it is available.  The solver is tried at tolerance 1e-3
   !       first; on failure it is retried at 1e-2 and then 3e-2 (ndgs stays
   !       2 -- raising it does not help, see the note below).  A RETRIED
   !       result is accepted only if it additionally reproduces the published
   !       Q_abs and Q_sca of this (a, lambda) to 1%.  The published tables
   !       come from an extended-precision run, so a "converged" solution of
   !       an ill-conditioned system that misses them by more than 1% is a
   !       false convergence and is refused.
   !
   !   regime 3, everything else (x > 50, or no tolerance converged)
   !       g of the VOLUME-EQUIVALENT SPHERE by Mie theory at the same m.
   !       This is an APPROXIMATION.  Its domain of validity is not a
   !       citation: the deviation of the sphere from the spheroid is measured
   !       on this very grid, on the converged points adjacent to the
   !       boundary, and the measured bound is printed by this program and
   !       written into the header of each G file.  It is the same class of
   !       stand-in Guillet et al. themselves used in the ultraviolet, where
   !       they interpolated to geometric optics using the volume-equivalent
   !       sphere's Mie albedo (their Eqs. 11-14, following Min et al. 2003).
   !
   ! WHY NOT spheroid_q.  That routine also selects a regime, but both of its
   ! asymptotic limits return g = 0.  For the Rayleigh limit that is correct
   ! physics; for the geometric-optics limit it is a placeholder its own
   ! source marks as unfit for a size distribution carrying large-grain
   ! weight, which is exactly this model's case.  This program therefore calls
   ! tmatrix_eval directly, inspects the verdict itself, and never takes a
   ! number from a failed solve.
   !
   ! WHY THE SOLVER FAILS.  tmatrix_core reports IERR = 3 when the multipole
   ! order reaches NPN1 = 100 without the convergence test being met.  At
   ! x ~ 6 the physical truncation order is ~15, so this is not a storage
   ! limit but the ill-conditioning of the extended-boundary-condition method
   ! for elongated particles in double precision: successive iterates never
   ! settle to within the tolerance and the order runs away.  Guillet et al.
   ! ran the extended-precision variant amplq.lp.f and reached a_V/lambda ~ 3
   ! for eps = 1/3 (their footnote 4); in double precision the boundary comes
   ! earlier, and where it comes is measured here rather than assumed.
   !
   ! COVERAGE OF THE OPTICAL CONSTANTS.  The BE index table ends at
   ! 0.0503 um, while LAMBDA.DAT starts at 0.04 um.  m is held at its
   ! 0.0503 um value over that gap; DustEM's own tables are extrapolated
   ! there too, and for this model only the _euv view of the product reaches
   ! those rows.  At the long-wavelength end nothing is extrapolated: every
   ! size on these grids has x < 0.1 beyond ~126 um, where regime 1 gives
   ! g = 0 without reading m at all.
   !
   ! Usage:
   !   ./spheroid_asymmetry_table.x            ! both populations
   !   ./spheroid_asymmetry_table.x amc        ! amCBE_0.3333x only
   !   ./spheroid_asymmetry_table.x sil        ! aSil2001BE6pctG_0.4x only
   !   ./spheroid_asymmetry_table.x amc test   ! every 8th size, for a smoke run
   !
   ! Written (relative to tmatrix/, where the Makefile drops the executable):
   !   ../data/g18d/oprop/G_<gtype>.DAT     the shipped table
   !   output/g18d_<gtype>_regimes.dat      diagnostics, not shipped: one row
   !                                        per (lambda, a) with x, a_V/lambda,
   !                                        the spheroid Q and g, the sphere g,
   !                                        the published Q, the regime and the
   !                                        tolerance that produced it.
   !
   ! `test` mode writes neither: it only prints the summary.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use tmatrix_api, only: tmatrix_options_t, tmatrix_result_t, &
                          tmatrix_workspace_t, tmatrix_workspace_init, &
                          tmatrix_workspace_finalize, tmatrix_eval
   use spheroid_optics, only: check_physical_bounds
   use read_index, only: refractive_index_t, load_index_energy_table, &
                         load_index_wavelength_table, refractive_index_at
   use effective_medium, only: maxwell_garnett_m
   use dustem_optical_tables, only: read_dustem_wavelength_grid, &
                                    read_dustem_sized_table
   use mie_mod, only: mie
   !$ use omp_lib
   implicit none

   integer, parameter :: wp = real64
   real(wp), parameter :: PI = acos(-1.0_wp)

   !! Upper bound of the Rayleigh dipole regime, the same 0.1 spheroid_optics
   !! uses: the neglected multipoles enter at O(x^2).
   real(wp), parameter :: X_RAYLEIGH_MAX = 0.1_wp
   !! Above this the T-matrix is not attempted at all.  It matches
   !! spheroid_optics::X_GEOMETRIC_MIN, and in practice the double-precision
   !! solver has given up an order of magnitude below it.
   real(wp), parameter :: X_TMATRIX_MAX = 50.0_wp
   !! Convergence tolerances tried in order.  1e-3 is the library default.
   integer,  parameter :: NTOL = 3
   real(wp), parameter :: TOLERANCE(NTOL) = [1.0e-3_wp, 1.0e-2_wp, 3.0e-2_wp]
   !! A retried solve must reproduce the published Q_abs and Q_sca this well.
   real(wp), parameter :: Q_REFERENCE_TOL = 0.01_wp
   !! Below this the sphere-vs-spheroid gap is not worth quoting as a bound on
   !! anything: regime 3 is never entered there.  The band that IS quoted is
   !! found from the table itself, as x >= the smallest x any regime-3 entry
   !! sits at.
   real(wp), parameter :: X_GAP_FLOOR = 0.1_wp

   integer, parameter :: REGIME_RAYLEIGH = 1, REGIME_TMATRIX = 2, REGIME_SPHERE = 3

   character(len=*), parameter :: DATA_DIR   = '../data/g18d/'
   !! Refractive-index tables, in the tree so that these G files can be
   !! regenerated from the release alone.  index_amcBE_ZMCB96 is the BE
   !! amorphous carbon of Zubko et al. (1996) in Draine's energy-ordered
   !! layout; eps_suvSil is his "smoothed UV" astrosilicate in the
   !! wavelength-ordered layout.  Both are copies of the files under
   !! Grain/opt/, recorded in data/g18d/where.txt.
   character(len=*), parameter :: DIEL_DIR   = '../data/dielectric/'
   character(len=*), parameter :: BE_INDEX_FILE     = DIEL_DIR//'index_amcBE_ZMCB96'
   character(len=*), parameter :: SUVSIL_INDEX_FILE = DIEL_DIR//'eps_suvSil'

   type :: population_t
      character(len=32)  :: tag         = ''
      character(len=64)  :: gtype       = ''
      real(wp)           :: eps_ba      = 1.0_wp
      logical            :: composite   = .false.  ! matrix + inclusions
      real(wp)           :: f_inclusion = 0.0_wp    ! volume fraction of inclusions
      real(wp)           :: a_model_max = 0.0_wp    ! [um], from GRAIN_G17_ModelD.DAT
      character(len=96)  :: material    = ''
   end type population_t

   type(population_t), parameter :: POPULATION(2) = [ &
      population_t('amc', 'amCBE_0.3333x', 1.0_wp/3.0_wp, .false., 0.0_wp, 1.540175_wp, &
                   'BE amorphous carbon (Zubko et al. 1996)'), &
      population_t('sil', 'aSil2001BE6pctG_0.4x', 0.4_wp, .true., 0.06_wp, 0.2918198_wp, &
                   'WD01 astrosilicate + 6% BE a-C, Maxwell Garnett') ]

   type(refractive_index_t) :: index_be, index_suvsil

   character(len=32) :: arg1, arg2
   logical :: test_mode
   integer :: ip, size_stride

   call load_index_energy_table(index_be, BE_INDEX_FILE)
   call load_index_wavelength_table(index_suvsil, SUVSIL_INDEX_FILE)

   arg1 = ''
   arg2 = ''
   if (command_argument_count() >= 1) call get_command_argument(1, arg1)
   if (command_argument_count() >= 2) call get_command_argument(2, arg2)
   test_mode = (trim(arg1) == 'test' .or. trim(arg2) == 'test')
   size_stride = 1
   if (test_mode) size_stride = 8

   do ip = 1, size(POPULATION)
      if (trim(arg1) /= '' .and. trim(arg1) /= 'test' .and. &
          trim(arg1) /= trim(POPULATION(ip)%tag)) cycle
      call build_one_population(POPULATION(ip), test_mode, size_stride)
   end do

contains

   subroutine build_one_population(pop, testing, stride)
      type(population_t), intent(in) :: pop
      logical, intent(in) :: testing
      integer, intent(in) :: stride

      real(wp), allocatable :: lambda(:), a_um(:), qblock(:,:,:)
      real(wp), allocatable :: nr(:), ki(:)
      logical,  allocatable :: m_known(:)
      real(wp), allocatable :: gtab(:,:), gsphere(:,:), qabs(:,:), qsca(:,:)
      integer,  allocatable :: regime(:,:), itol(:,:)
      integer  :: nwave, nsize, nblock, iw, n_clamped
      real(wp) :: lam_min_be, lam_min_mix
      real(wp) :: t_start, t_end
      character(len=512) :: qpath

      write(*,'(a)') repeat('-', 78)
      write(*,'(a,a)') ' population: ', trim(pop%gtype)

      call read_dustem_wavelength_grid(DATA_DIR//'oprop/LAMBDA.DAT', nwave, lambda)
      qpath = DATA_DIR//'oprop/Q_'//trim(pop%gtype)//'.DAT'
      nblock = 2
      call read_dustem_sized_table(trim(qpath), nwave, nsize, a_um, qblock, nblock)
      if (nblock /= 2) then
         write(error_unit,'(a,a)') ' expected two Q blocks in ', trim(qpath)
         stop 1
      end if
      write(*,'(a,i0,a,i0,a,es10.3,a,es10.3,a)') '   grid: ', nwave, ' wavelengths x ', &
           nsize, ' radii (', a_um(1), ' to ', a_um(nsize), ' um)'

      ! ---- refractive index on the wavelength grid --------------------------
      ! m is needed only where some radius reaches x >= X_RAYLEIGH_MAX, i.e.
      ! lambda <= 2 pi a_max / X_RAYLEIGH_MAX.  Beyond that every entry is
      ! regime 1 and the index is never consulted, which is what keeps the
      ! tables' long-wavelength ends (1985 um for BE, 1000 um for suvSil) from
      ! ever being extrapolated.
      allocate(nr(nwave), ki(nwave), m_known(nwave))
      nr = 0.0_wp;  ki = 0.0_wp;  m_known = .false.
      lam_min_be  = index_be%lambda(1)
      lam_min_mix = max(index_be%lambda(1), index_suvsil%lambda(1))
      n_clamped = 0
      do iw = 1, nwave
         if (2.0_wp*PI*a_um(nsize)/lambda(iw) < X_RAYLEIGH_MAX) cycle
         call refractive_index_of(pop, lambda(iw), nr(iw), ki(iw), lam_min_be, &
                                  lam_min_mix, n_clamped)
         m_known(iw) = .true.
      end do
      if (n_clamped > 0) write(*,'(a,i0,a,es10.3,a)') &
         '   m held constant over ', n_clamped, &
         ' wavelengths shortward of the index table (', &
         merge(lam_min_be, lam_min_mix, .not. pop%composite), ' um)'

      ! ---- the sweep --------------------------------------------------------
      allocate(gtab(nwave, nsize), gsphere(nwave, nsize))
      allocate(qabs(nwave, nsize), qsca(nwave, nsize))
      allocate(regime(nwave, nsize), itol(nwave, nsize))
      gtab = 0.0_wp;  gsphere = 0.0_wp
      qabs = 0.0_wp;  qsca = 0.0_wp
      regime = REGIME_RAYLEIGH;  itol = 0

      call wall_clock(t_start)
      !$omp parallel default(shared)
      call sweep_sizes(pop, lambda, a_um, nr, ki, m_known, qblock, stride, &
                       gtab, gsphere, qabs, qsca, regime, itol)
      !$omp end parallel
      call wall_clock(t_end)

      write(*,'(a,f9.1,a)') '   wall clock: ', t_end - t_start, ' s'
      !$omp parallel
      !$omp masked
      !$ write(*,'(a,i0)') '   OpenMP threads: ', omp_get_num_threads()
      !$omp end masked
      !$omp end parallel

      call report_population(pop, lambda, a_um, gtab, gsphere, qabs, qsca, &
                             qblock, regime, itol, stride)

      if (.not. testing) then
         call write_gfile(pop, lambda, a_um, gtab, gsphere, regime, trim(qpath))
         call write_diagnostics(pop, lambda, a_um, gtab, gsphere, qabs, qsca, &
                                qblock, regime, itol)
         call write_convergence_boundary(pop, lambda, a_um, regime, itol)
      else
         write(*,'(a)') '   test mode: no file written'
      end if

      ! Sizes not visited in test mode keep g = 0; say so rather than letting
      ! the summary above look like a full table.
      if (stride > 1) write(*,'(a,i0,a)') &
         '   test mode: every ', stride, 'th radius was evaluated'

      deallocate(lambda, a_um, qblock, nr, ki, m_known)
      deallocate(gtab, gsphere, qabs, qsca, regime, itol)
   end subroutine build_one_population


   subroutine sweep_sizes(pop, lambda, a_um, nr, ki, m_known, qblock, stride, &
                          gtab, gsphere, qabs, qsca, regime, itol)
      !! One OpenMP thread owns one T-matrix workspace and takes whole radii.
      type(population_t), intent(in) :: pop
      real(wp), intent(in) :: lambda(:), a_um(:), nr(:), ki(:)
      logical,  intent(in) :: m_known(:)
      real(wp), intent(in) :: qblock(:,:,:)
      integer,  intent(in) :: stride
      real(wp), intent(inout) :: gtab(:,:), gsphere(:,:), qabs(:,:), qsca(:,:)
      integer,  intent(inout) :: regime(:,:), itol(:,:)

      type(tmatrix_workspace_t) :: work
      type(tmatrix_options_t)   :: opt
      type(tmatrix_result_t)    :: res
      integer  :: ia, iw, k, nwave, nsize
      real(wp) :: x, qe, qs, qa, alb, g
      logical  :: ok, accepted
      character(len=96) :: reason

      nwave = size(lambda)
      nsize = size(a_um)
      call tmatrix_workspace_init(work)
      opt%shape        = -1
      opt%ndgs         = 2
      opt%aspect_ratio = pop%eps_ba

      !$omp do schedule(dynamic,1)
      do ia = 1, nsize
         if (mod(ia-1, stride) /= 0) cycle
         do iw = 1, nwave
            x = 2.0_wp*PI*a_um(ia)/lambda(iw)
            if (x < X_RAYLEIGH_MAX) then
               regime(iw, ia) = REGIME_RAYLEIGH
               gtab(iw, ia)   = 0.0_wp
               cycle
            end if
            if (.not. m_known(iw)) then
               write(error_unit,'(a)') ' sweep_sizes: m needed where it was not evaluated.'
               stop 1
            end if

            ! Volume-equivalent sphere, always: it is the regime-3 answer and
            ! the measurement of how far the sphere sits from the spheroid.
            call mie(nr(iw), ki(iw), x, qe, qs, qa, alb, g)
            gsphere(iw, ia) = g

            accepted = .false.
            if (x <= X_TMATRIX_MAX) then
               do k = 1, NTOL
                  opt%tolerance = TOLERANCE(k)
                  call tmatrix_eval(work, a_um(ia), lambda(iw), nr(iw), ki(iw), opt, res)
                  if (res%status /= 0 .or. res%legacy_ierr /= 0) cycle
                  call check_physical_bounds(res%qext, res%qabs, res%qsca, &
                                             res%albedo, res%asymmetry, ok, reason)
                  if (.not. ok) cycle
                  if (k > 1) then
                     if (.not. matches_published_q(res%qabs, res%qsca, &
                              qblock(iw, ia, 1), qblock(iw, ia, 2))) cycle
                  end if
                  regime(iw, ia) = REGIME_TMATRIX
                  itol(iw, ia)   = k
                  gtab(iw, ia)   = res%asymmetry
                  qabs(iw, ia)   = res%qabs
                  qsca(iw, ia)   = res%qsca
                  accepted = .true.
                  exit
               end do
            end if

            if (.not. accepted) then
               regime(iw, ia) = REGIME_SPHERE
               itol(iw, ia)   = 0
               gtab(iw, ia)   = gsphere(iw, ia)
               qabs(iw, ia)   = qa
               qsca(iw, ia)   = qs
            end if
         end do
      end do
      !$omp end do
      call tmatrix_workspace_finalize(work)
   end subroutine sweep_sizes


   logical function matches_published_q(qa, qs, qa_ref, qs_ref) result(ok)
      !! The published tables come from an extended-precision solve.  A
      !! retried, looser-tolerance result is trusted only where it reproduces
      !! them.  A non-positive reference cannot decide anything, so the test
      !! fails and the entry falls to the sphere stand-in.
      real(wp), intent(in) :: qa, qs, qa_ref, qs_ref
      ok = .false.
      if (qa_ref <= 0.0_wp .or. qs_ref <= 0.0_wp) return
      if (abs(qa - qa_ref)/qa_ref > Q_REFERENCE_TOL) return
      if (abs(qs - qs_ref)/qs_ref > Q_REFERENCE_TOL) return
      ok = .true.
   end function matches_published_q


   subroutine refractive_index_of(pop, lam, nr_out, ki_out, lam_min_be, &
                                  lam_min_mix, n_clamped)
      !! m of the population's material at one wavelength, with the composite
      !! mixed here rather than in the table.
      type(population_t), intent(in) :: pop
      real(wp), intent(in)  :: lam, lam_min_be, lam_min_mix
      real(wp), intent(out) :: nr_out, ki_out
      integer,  intent(inout) :: n_clamped
      real(wp) :: lam_use, nr_m, ki_m, nr_i, ki_i

      if (.not. pop%composite) then
         lam_use = max(lam, lam_min_be)
         if (lam < lam_min_be) n_clamped = n_clamped + 1
         call refractive_index_at(index_be, lam_use, nr_out, ki_out)
      else
         lam_use = max(lam, lam_min_mix)
         if (lam < lam_min_mix) n_clamped = n_clamped + 1
         call refractive_index_at(index_suvsil, lam_use, nr_m, ki_m)
         call refractive_index_at(index_be,     lam_use, nr_i, ki_i)
         call maxwell_garnett_m(nr_m, ki_m, nr_i, ki_i, pop%f_inclusion, nr_out, ki_out)
      end if
   end subroutine refractive_index_of


   subroutine report_population(pop, lambda, a_um, gtab, gsphere, qabs, qsca, &
                                qblock, regime, itol, stride)
      !! Everything the report needs: regime census, the convergence boundary
      !! in x and a_V/lambda, the agreement with the published Q where the
      !! T-matrix converged, and the sphere-vs-spheroid bound on g.
      type(population_t), intent(in) :: pop
      real(wp), intent(in) :: lambda(:), a_um(:), gtab(:,:), gsphere(:,:)
      real(wp), intent(in) :: qabs(:,:), qsca(:,:), qblock(:,:,:)
      integer,  intent(in) :: regime(:,:), itol(:,:), stride

      integer  :: iw, ia, k, nwave, nsize, nvis, nmod
      integer  :: n_reg(3), n_reg_mod(3), n_tol(NTOL), n_gneg
      real(wp) :: xmax_tol(NTOL), avlam_max(NTOL)
      real(wp) :: x, dq_a, dq_s, dqa_max, dqs_max
      real(wp) :: sum_a, sum_s, gmin, gmax
      integer  :: n_dq

      nwave = size(lambda);  nsize = size(a_um)
      n_reg = 0;  n_reg_mod = 0;  n_tol = 0;  n_gneg = 0
      xmax_tol = 0.0_wp;  avlam_max = 0.0_wp
      dqa_max = 0.0_wp;  dqs_max = 0.0_wp;  sum_a = 0.0_wp;  sum_s = 0.0_wp
      n_dq = 0
      gmin = huge(1.0_wp);  gmax = -huge(1.0_wp)
      nvis = 0;  nmod = 0

      do ia = 1, nsize
         if (mod(ia-1, stride) /= 0) cycle
         do iw = 1, nwave
            nvis = nvis + 1
            k = regime(iw, ia)
            n_reg(k) = n_reg(k) + 1
            if (a_um(ia) <= pop%a_model_max) then
               nmod = nmod + 1
               n_reg_mod(k) = n_reg_mod(k) + 1
            end if
            gmin = min(gmin, gtab(iw, ia));  gmax = max(gmax, gtab(iw, ia))
            if (gtab(iw, ia) < 0.0_wp) n_gneg = n_gneg + 1
            x = 2.0_wp*PI*a_um(ia)/lambda(iw)
            if (k == REGIME_TMATRIX) then
               n_tol(itol(iw,ia)) = n_tol(itol(iw,ia)) + 1
               xmax_tol(itol(iw,ia))  = max(xmax_tol(itol(iw,ia)), x)
               avlam_max(itol(iw,ia)) = max(avlam_max(itol(iw,ia)), a_um(ia)/lambda(iw))
               ! Agreement with the published table, and the sphere gap, both
               ! measured only where the T-matrix produced the answer.
               if (qblock(iw,ia,1) > 0.0_wp .and. qblock(iw,ia,2) > 0.0_wp) then
                  dq_a = abs(qabs(iw,ia) - qblock(iw,ia,1))/qblock(iw,ia,1)
                  dq_s = abs(qsca(iw,ia) - qblock(iw,ia,2))/qblock(iw,ia,2)
                  dqa_max = max(dqa_max, dq_a);  dqs_max = max(dqs_max, dq_s)
                  sum_a = sum_a + dq_a;  sum_s = sum_s + dq_s;  n_dq = n_dq + 1
               end if
            end if
         end do
      end do

      write(*,'(a)') '   regime census (whole table / within the model size range):'
      write(*,'(a,i8,a,f6.2,a,i8,a,f6.2,a)') &
         '     1 Rayleigh g=0     ', n_reg(1), ' (', 100.0_wp*n_reg(1)/max(nvis,1), '%)  ', &
         n_reg_mod(1), ' (', 100.0_wp*n_reg_mod(1)/max(nmod,1), '%)'
      write(*,'(a,i8,a,f6.2,a,i8,a,f6.2,a)') &
         '     2 T-matrix         ', n_reg(2), ' (', 100.0_wp*n_reg(2)/max(nvis,1), '%)  ', &
         n_reg_mod(2), ' (', 100.0_wp*n_reg_mod(2)/max(nmod,1), '%)'
      write(*,'(a,i8,a,f6.2,a,i8,a,f6.2,a)') &
         '     3 sphere stand-in  ', n_reg(3), ' (', 100.0_wp*n_reg(3)/max(nvis,1), '%)  ', &
         n_reg_mod(3), ' (', 100.0_wp*n_reg_mod(3)/max(nmod,1), '%)'
      do k = 1, NTOL
         write(*,'(a,es9.2,a,i8,a,f8.3,a,f8.3)') &
            '     tolerance ', TOLERANCE(k), ': ', n_tol(k), &
            ' accepted, x_max = ', xmax_tol(k), ', (a_V/lambda)_max = ', avlam_max(k)
      end do
      if (n_dq > 0) then
         write(*,'(a,f8.4,a,f8.4,a)') '   vs published Q where the T-matrix answered: ' // &
            'max |dQabs| = ', 100.0_wp*dqa_max, '%, mean ', 100.0_wp*sum_a/n_dq, '%'
         write(*,'(a,f8.4,a,f8.4,a)') '                                              ' // &
            'max |dQsca| = ', 100.0_wp*dqs_max, '%, mean ', 100.0_wp*sum_s/n_dq, '%'
      end if
      call report_sphere_gap_bins(pop, lambda, a_um, gtab, gsphere, regime, stride)
      write(*,'(a,f9.5,a,f9.5,a,i0)') '   g range [', gmin, ', ', gmax, &
           '],  entries with g < 0: ', n_gneg
   end subroutine report_population



   subroutine report_sphere_gap_bins(pop, lambda, a_um, gtab, gsphere, regime, stride)
      !! How far the volume-equivalent sphere sits from the spheroid in g, as a
      !! function of x, over the points where the T-matrix gave the spheroid
      !! answer.  Binned because the two are not uniformly close: they agree to
      !! a few parts in a thousand once x is large, and part company most
      !! around x ~ 1.  Regime 3 is only ever entered at large x, so the last
      !! bins are the ones that bound the approximation.
      type(population_t), intent(in) :: pop
      real(wp), intent(in) :: lambda(:), a_um(:), gtab(:,:), gsphere(:,:)
      integer,  intent(in) :: regime(:,:), stride
      integer,  parameter  :: NBIN = 6
      real(wp), parameter  :: EDGE(NBIN+1) = &
         [0.1_wp, 0.5_wp, 1.0_wp, 2.0_wp, 5.0_wp, 10.0_wp, huge(1.0_wp)]
      integer  :: cnt(NBIN), cnt_m(NBIN), ia, iw, ib
      real(wp) :: dmax(NBIN), dsum(NBIN), dmax_m(NBIN), x, dg

      cnt = 0;  dmax = 0.0_wp;  dsum = 0.0_wp
      cnt_m = 0;  dmax_m = 0.0_wp
      do ia = 1, size(a_um)
         if (mod(ia-1, stride) /= 0) cycle
         do iw = 1, size(lambda)
            if (regime(iw,ia) /= REGIME_TMATRIX) cycle
            x  = 2.0_wp*PI*a_um(ia)/lambda(iw)
            dg = abs(gsphere(iw,ia) - gtab(iw,ia))
            do ib = 1, NBIN
               if (x >= EDGE(ib) .and. x < EDGE(ib+1)) then
                  cnt(ib)  = cnt(ib) + 1
                  dsum(ib) = dsum(ib) + dg
                  dmax(ib) = max(dmax(ib), dg)
                  if (a_um(ia) <= pop%a_model_max) then
                     cnt_m(ib)  = cnt_m(ib) + 1
                     dmax_m(ib) = max(dmax_m(ib), dg)
                  end if
                  exit
               end if
            end do
         end do
      end do
      write(*,'(a)') '   |g(sphere) - g(spheroid)| over regime 2, by size parameter'
      write(*,'(a)') '   (the second column is restricted to the radii the model reaches):'
      do ib = 1, NBIN
         if (cnt(ib) == 0) cycle
         if (ib < NBIN) then
            write(*,'(a,f6.2,a,f6.2,a,i7,a,f8.5,a,f8.5,a,i7,a,f8.5)') '     ', EDGE(ib), &
                 ' <= x < ', EDGE(ib+1), ':  n = ', cnt(ib), ',  mean ', dsum(ib)/cnt(ib), &
                 ',  max ', dmax(ib), '   | in-model n = ', cnt_m(ib), ',  max ', dmax_m(ib)
         else
            write(*,'(a,f6.2,a,i7,a,f8.5,a,f8.5,a,i7,a,f8.5)') '     ', EDGE(ib), &
                 ' <= x       :  n = ', cnt(ib), ',  mean ', dsum(ib)/cnt(ib), &
                 ',  max ', dmax(ib), '   | in-model n = ', cnt_m(ib), ',  max ', dmax_m(ib)
         end if
      end do
   end subroutine report_sphere_gap_bins


   subroutine write_convergence_boundary(pop, lambda, a_um, regime, itol)
      !! The T-matrix convergence boundary of this population, size by size:
      !! the largest x (shortest lambda) at which a solve was accepted, and the
      !! tolerance it took.  This is what Guillet et al. footnote 4 quotes as
      !! a_V/lambda for their extended-precision build.
      type(population_t), intent(in) :: pop
      real(wp), intent(in) :: lambda(:), a_um(:)
      integer,  intent(in) :: regime(:,:), itol(:,:)
      integer  :: u, ia, iw, iw_best, k_best
      real(wp) :: x, x_best
      character(len=512) :: path

      path = 'output/g18d_'//trim(pop%gtype)//'_boundary.dat'
      open(newunit=u, file=trim(path), status='replace', action='write')
      write(u,'(a,a)') '# T-matrix convergence boundary for ', trim(pop%gtype)
      write(u,'(a)') '# Written by tmatrix/driver/spheroid_asymmetry_table.f90.  For each'
      write(u,'(a)') '# radius: the largest size parameter at which a solve was accepted, the'
      write(u,'(a)') '# wavelength and a_V/lambda there, and the tolerance it took.  A radius'
      write(u,'(a)') '# with no accepted solve carries zeros.'
      write(u,'(a)') '#     a_V[um]       x_max  lambda[um]     aV/lam    tol      in_model'
      do ia = 1, size(a_um)
         x_best = 0.0_wp;  iw_best = 0;  k_best = 0
         do iw = 1, size(lambda)
            if (regime(iw,ia) /= REGIME_TMATRIX) cycle
            x = 2.0_wp*PI*a_um(ia)/lambda(iw)
            if (x > x_best) then
               x_best = x;  iw_best = iw;  k_best = itol(iw,ia)
            end if
         end do
         if (iw_best > 0) then
            write(u,'(4es12.4,es10.2,i8)') a_um(ia), x_best, lambda(iw_best), &
                 a_um(ia)/lambda(iw_best), TOLERANCE(k_best), &
                 merge(1, 0, a_um(ia) <= pop%a_model_max)
         else
            write(u,'(4es12.4,es10.2,i8)') a_um(ia), 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                 merge(1, 0, a_um(ia) <= pop%a_model_max)
         end if
      end do
      close(u)
      write(*,'(a,a)') '   wrote ', trim(path)
   end subroutine write_convergence_boundary

   subroutine write_gfile(pop, lambda, a_um, gtab, gsphere, regime, qpath)
      type(population_t), intent(in) :: pop
      real(wp), intent(in) :: lambda(:), a_um(:), gtab(:,:), gsphere(:,:)
      integer,  intent(in) :: regime(:,:)
      character(len=*), intent(in) :: qpath
      integer :: u, iw, ia, i0, i1, nwave, nsize, n_reg(3)
      real(wp) :: dg_max, x_at, x3_min, dg_band, x_band, dg_model
      character(len=512) :: path
      character(len=8)  :: date
      character(len=10) :: time

      nwave = size(lambda);  nsize = size(a_um)
      n_reg = 0
      do ia = 1, nsize
         do iw = 1, nwave
            n_reg(regime(iw,ia)) = n_reg(regime(iw,ia)) + 1
         end do
      end do
      call sphere_gap(pop, lambda, a_um, gtab, gsphere, regime, &
                      dg_max, x_at, x3_min, dg_band, x_band, dg_model)
      call date_and_time(date=date, time=time)

      path = DATA_DIR//'oprop/G_'//trim(pop%gtype)//'.DAT'
      open(newunit=u, file=trim(path), status='replace', action='write')
      write(u,'(a)') '# g=<cos(theta)>  scattering factors for DUSTEM'
      write(u,'(a)') '#'
      write(u,'(a)') '# NOT FROM THE DustEM DISTRIBUTION.  The distribution ships no G_ file'
      write(u,'(a)') '# for this population; Guillet et al. (2018, A&A 610, A16), whose model D'
      write(u,'(a)') '# this is, computed no asymmetry parameter, and DustEM reads a G_ file'
      write(u,'(a)') '# only under its `pdr` run keyword, which this model does not set.  This'
      write(u,'(a)') '# table was computed by SEDust with'
      write(u,'(a)') '#   tmatrix/driver/spheroid_asymmetry_table.f90'
      write(u,'(a,a,a,a)') '# on ', date(1:4)//'-'//date(5:6)//'-'//date(7:8), ' ', &
                            time(1:2)//':'//time(3:4)//' local time.'
      write(u,'(a)') '#'
      write(u,'(a,a)') '# material : ', trim(pop%material)
      if (pop%composite) then
         write(u,'(a,a)') '#             matrix    ', SUVSIL_INDEX_FILE
         write(u,'(a,a)') '#             inclusion ', BE_INDEX_FILE
         write(u,'(a,f6.3,a)') '#             volume fraction of inclusions ', &
              pop%f_inclusion, ', Maxwell Garnett'
         write(u,'(a)') '#             (Bohren & Huffman 1983 Sect. 8.5; Guillet et al. Eq. 28)'
      else
         write(u,'(a,a)') '#             ', BE_INDEX_FILE
      end if
      write(u,'(a,f8.5,a)') '# shape    : prolate spheroid, axis ratio b/a = ', pop%eps_ba, &
           ' (Guillet et al. Table 3, model D)'
      write(u,'(a)') '# geometry : random orientation.  Q and g are defined against the'
      write(u,'(a)') '#            volume-equivalent sphere radius a_V, as Guillet et al.'
      write(u,'(a)') '#            Sect. 3.6 defines them, and the radii below are a_V.'
      write(u,'(a,a)') '# radii    : the same grid as ', trim(qpath)
      write(u,'(a)') '# rows     : the wavelengths of oprop/LAMBDA.DAT, in that order.'
      write(u,'(a)') '#'
      write(u,'(a)') '# METHOD, by size parameter x = 2 pi a_V / lambda:'
      write(u,'(a)') '#   x < 0.1              g = 0.  Rayleigh dipole scattering is'
      write(u,'(a)') '#                        forward-backward symmetric; the correction'
      write(u,'(a)') '#                        is O(x^2).  Exact, not an approximation.'
      write(u,'(a)') '#   0.1 <= x <= 50,      g = the random-orientation T-matrix asymmetry'
      write(u,'(a)') '#   solver converged     (Mishchenko), which is the analytic form of'
      write(u,'(a)') '#                        the 768-direction HEALPix orientation average'
      write(u,'(a)') '#                        Guillet et al. Eq. 19 takes.  Tolerance 1e-3,'
      write(u,'(a)') '#                        retried at 1e-2 and 3e-2; a retried solve is'
      write(u,'(a)') '#                        kept only where it also reproduces the'
      write(u,'(a)') '#                        published Q_abs and Q_sca to 1%.'
      write(u,'(a)') '#   otherwise            g of the VOLUME-EQUIVALENT SPHERE by Mie'
      write(u,'(a)') '#                        theory at the same m.  APPROXIMATION.  In'
      write(u,'(a)') '#                        double precision the T-matrix gives up well'
      write(u,'(a)') '#                        below the x = 50 bound (Guillet et al. ran an'
      write(u,'(a)') '#                        extended-precision build and reached'
      write(u,'(a)') '#                        a_V/lambda ~ 3 for b/a = 1/3, their'
      write(u,'(a)') '#                        footnote 4), so this covers the large-x end.'
      write(u,'(a)') '#                        MEASURED bound, from the regime-2 cells of'
      write(u,'(a)') '#                        this very table, where both forms are known:'
      write(u,'(a,f7.2,a)') '#                          the stand-in is entered only at x >= ', &
           x3_min, ', and over'
      write(u,'(a,f7.4,a,f7.2,a)') '#                          that band |g_sphere - g_spheroid| <= ', &
           dg_band, ' (at x = ', x_band, '),'
      write(u,'(a,f7.4,a)') '#                          and <= ', dg_model, &
           ' over the radii the model reaches.'
      write(u,'(a,f7.4,a,f7.2,a)') '#                        Anywhere in regime 2 the gap reaches ', &
           dg_max, ' (at x = ', x_at, '),'
      write(u,'(a)') '#                        but that x is below the threshold above, so the'
      write(u,'(a)') '#                        stand-in is never reached there: the two forms'
      write(u,'(a)') '#                        part company at moderate x and converge as x'
      write(u,'(a)') '#                        grows.'
      write(u,'(a,i0,a,i0,a,i0,a)') '# census  : ', n_reg(1), ' Rayleigh, ', n_reg(2), &
           ' T-matrix, ', n_reg(3), ' sphere entries.'
      if (.not. pop%composite) then
         write(u,'(a)') '# note    : the BE index table ends at 0.0503 um; m is held constant'
         write(u,'(a)') '#           from there to the 0.04 um end of LAMBDA.DAT.  The Q_ table'
         write(u,'(a)') '#           of this population is identically zero for lambda < 0.09 um,'
         write(u,'(a)') '#           so those rows carry no weight in this model.'
      end if
      write(u,'(a)') '#'
      write(u,'(a)') '# nsize (number of sampled sizes in this file)'
      write(u,'(a)') '# sizes (microns)'
      write(u,'(a)') '# g values: g(size1),...., g(sizeN)'
      write(u,'(a)') '#'
      write(u,'(i12)') nsize
      do i0 = 1, nsize, 10
         i1 = min(i0+9, nsize)
         write(u,'(*(1x,es11.4))') a_um(i0:i1)
      end do
      write(u,'(a)') '#### g-factor ####              (lines are lambda from LAMBDA.DAT and columns are the N sizes)'
      do iw = 1, nwave
         write(u,'(*(1x,es11.4))') gtab(iw, 1:nsize)
      end do
      close(u)
      write(*,'(a,a)') '   wrote ', trim(path)
   end subroutine write_gfile


   subroutine sphere_gap(pop, lambda, a_um, gtab, gsphere, regime, &
                         dg_max, x_at, x3_min, dg_band, x_band, dg_model)
      !! |g(sphere) - g(spheroid)| over the points where the T-matrix produced
      !! the spheroid answer -- the whole measurable extent of the regime-3
      !! approximation, since outside regime 2 there is nothing to compare the
      !! sphere against.
      !!
      !! Three numbers come back.  dg_max is the largest gap anywhere in
      !! regime 2; it falls around x ~ 1, where regime 3 is never entered, so
      !! it bounds nothing this table does.  dg_band is the largest gap over
      !! x >= x3_min, where x3_min is the smallest size parameter any regime-3
      !! entry sits at: that IS the band the stand-in operates in, and dg_band
      !! is the bound on it.  dg_model is dg_band restricted to the radii the
      !! model's size distribution reaches.
      type(population_t), intent(in) :: pop
      real(wp), intent(in) :: lambda(:), a_um(:), gtab(:,:), gsphere(:,:)
      integer,  intent(in) :: regime(:,:)
      real(wp), intent(out) :: dg_max, x_at, x3_min, dg_band, x_band, dg_model
      integer :: iw, ia
      real(wp) :: dg, x

      x3_min = huge(1.0_wp)
      do ia = 1, size(a_um)
         do iw = 1, size(lambda)
            if (regime(iw,ia) /= REGIME_SPHERE) cycle
            x3_min = min(x3_min, 2.0_wp*PI*a_um(ia)/lambda(iw))
         end do
      end do
      if (x3_min > huge(1.0_wp)*0.5_wp) x3_min = X_GAP_FLOOR

      dg_max = 0.0_wp;  x_at = 0.0_wp
      dg_band = 0.0_wp; x_band = 0.0_wp;  dg_model = 0.0_wp
      do ia = 1, size(a_um)
         do iw = 1, size(lambda)
            if (regime(iw,ia) /= REGIME_TMATRIX) cycle
            x = 2.0_wp*PI*a_um(ia)/lambda(iw)
            dg = abs(gsphere(iw,ia) - gtab(iw,ia))
            if (dg > dg_max) then
               dg_max = dg;  x_at = x
            end if
            if (x < x3_min) cycle
            if (dg > dg_band) then
               dg_band = dg;  x_band = x
            end if
            if (a_um(ia) <= pop%a_model_max) dg_model = max(dg_model, dg)
         end do
      end do
   end subroutine sphere_gap


   subroutine write_diagnostics(pop, lambda, a_um, gtab, gsphere, qabs, qsca, &
                                qblock, regime, itol)
      type(population_t), intent(in) :: pop
      real(wp), intent(in) :: lambda(:), a_um(:), gtab(:,:), gsphere(:,:)
      real(wp), intent(in) :: qabs(:,:), qsca(:,:), qblock(:,:,:)
      integer,  intent(in) :: regime(:,:), itol(:,:)
      integer :: u, iw, ia
      real(wp) :: tolv
      character(len=512) :: path

      path = 'output/g18d_'//trim(pop%gtype)//'_regimes.dat'
      open(newunit=u, file=trim(path), status='replace', action='write')
      write(u,'(a,a)') '# regime diagnostics for G_', trim(pop%gtype)//'.DAT'
      write(u,'(a)') '# Written by tmatrix/driver/spheroid_asymmetry_table.f90.  NOT shipped;'
      write(u,'(a)') '# this is the record of how every entry of the G table was obtained.'
      write(u,'(a)') '# regime 1 = Rayleigh dipole (g = 0), 2 = random-orientation T-matrix,'
      write(u,'(a)') '# 3 = volume-equivalent sphere by Mie.  tol is the T-matrix convergence'
      write(u,'(a)') '# tolerance that was accepted (0 outside regime 2).  Qabs_pub / Qsca_pub'
      write(u,'(a)') '# are the DustEM published values at the same cell.'
      write(u,'(a)') '#  lambda[um]     a_V[um]           x      a_V/lam      Q_abs      ' // &
                     'Q_sca          g    g_sphere   Qabs_pub   Qsca_pub  reg        tol'
      do ia = 1, size(a_um)
         do iw = 1, size(lambda)
            tolv = 0.0_wp
            if (itol(iw,ia) > 0) tolv = TOLERANCE(itol(iw,ia))
            write(u,'(2es13.5,2es13.5,6es11.3,i5,es11.2)') &
               lambda(iw), a_um(ia), 2.0_wp*PI*a_um(ia)/lambda(iw), a_um(ia)/lambda(iw), &
               qabs(iw,ia), qsca(iw,ia), gtab(iw,ia), gsphere(iw,ia), &
               qblock(iw,ia,1), qblock(iw,ia,2), regime(iw,ia), tolv
         end do
      end do
      close(u)
      write(*,'(a,a)') '   wrote ', trim(path)
   end subroutine write_diagnostics


   subroutine wall_clock(t)
      real(wp), intent(out) :: t
      integer :: c, r
      call system_clock(count=c, count_rate=r)
      t = real(c, wp)/real(max(r,1), wp)
   end subroutine wall_clock

end program spheroid_asymmetry_table
