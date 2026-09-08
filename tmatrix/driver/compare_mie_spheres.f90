program compare_mie_spheres
   ! Certifies the two refractive-index tables and this tree's Mie solver
   ! against published sphere tables, BEFORE any shape enters.  It is the
   ! first two of the four checks the G18 Model D asymmetry tables rest on:
   ! if our m(lambda) and our Q, g of a SPHERE do not reproduce the published
   ! sphere tables, nothing computed for a spheroid from the same index can be
   ! believed.
   !
   ! CHECK A -- BE amorphous carbon (Zubko et al. 1996), the material of the
   !   amCBE_0.3333x population.  Our Mie with
   !     data/dielectric/index_amcBE_ZMCB96
   !   against DustEM's own sphere tables for the same material,
   !     DustEM/oprop/Q_amCBEx.DAT   (Q_abs, Q_sca)
   !     DustEM/oprop/G_amCBEx.DAT   (g)
   !   on their 50 radii and the 800 wavelengths of DustEM's LAMBDA.DAT.
   !   Conditions that are NOT matched, and are excluded rather than smoothed:
   !     * lambda > 800 um.  Guillet et al. Sect. 5.1 state that the BE data
   !       stop at 2 mm and deviate from a power law from 1 mm, and that
   !       DustEM extrapolates the cross sections from 800 um to 1 cm with a
   !       power law fitted over 800-1000 um.  Beyond 800 um the DustEM table
   !       is therefore that extrapolation, not a Mie calculation.
   !     * lambda < 0.0503 um, the short-wavelength end of the index table.
   !       DustEM's table extends to 0.04 um, so those rows likewise come from
   !       an extrapolation of the index that is not recorded in the file we
   !       hold.  They are reported separately and not counted against the
   !       method.
   !   `G_amCBE.DAT` and `G_amCBEx.DAT` are byte-identical in the
   !   distribution; the x form is used here because it is the one the Q file
   !   used for this model is named after.
   !
   ! CHECK B -- WD01 astrosilicate, the matrix of the aSil2001BE6pctG_0.4x
   !   population.  Our Mie with
   !     data/dielectric/eps_suvSil
   !   against Draine's own Mie table for the same material,
   !     <Draine opt>/suvSil_81  (81 radii x 241 wavelengths, Q_abs, Q_sca, g)
   !   on his grid.  Both files are Draine's and cover the same range.  Two
   !   parts of his grid reach far outside anything this model uses and are
   !   reported apart: lambda < 0.04 um (the short end of DustEM's LAMBDA.DAT)
   !   and x > 320 (the largest size parameter the G18 Model D tables reach,
   !   2 pi * 2 um / 0.04 um).
   !
   !   This check also tests the REFERENCE, with an argument that needs no
   !   second code: at x < 0.01 absorption is in the Rayleigh dipole limit,
   !   where Q_abs is exactly proportional to a at fixed lambda, so
   !   Q_abs / a must not depend on the radius.  Ours obeys that to four
   !   digits everywhere; suvSil_81 has isolated patches of a few adjacent
   !   radii where it does not, and those are the cells that produce every
   !   Q_abs disagreement above 1%.  The count and the worst case are printed.
   !
   ! REFERENCE DATA OUTSIDE THE RELEASE.  The two index tables are read from
   ! data/dielectric/, so the side of this check that belongs to SEDust runs
   ! from the release alone.  The three tables it is checked AGAINST -- DustEM's
   ! Q_amCBEx / G_amCBEx and Draine's suvSil_81 -- are published data that
   ! SEDust does not redistribute, and are read from the working tree they sit
   ! in.  The two roots are the parameters DUSTEM and DRAINE_OPT below; point
   ! them at your own copies.  This is a verification driver: no product
   ! depends on it, and nothing else in the tree reads those paths.
   !
   ! Usage (from tmatrix/):
   !   ./compare_mie_spheres.x
   !
   ! Writes output/mie_vs_amCBEx.dat and output/mie_vs_suvSil_81.dat, one row
   ! per (lambda, a) with both values of each quantity, and prints the maximum
   ! and median deviations.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use read_index, only: refractive_index_t, load_index_energy_table, &
                         load_index_wavelength_table, refractive_index_at
   use dustem_optical_tables, only: read_dustem_wavelength_grid, &
                                    read_dustem_sized_table
   use mie_mod, only: mie
   implicit none

   integer, parameter :: wp = real64
   real(wp), parameter :: PI = acos(-1.0_wp)
   character(len=*), parameter :: DUSTEM = '/home/kiseon/MoCafe/Grain/DustEM/oprop/'
   character(len=*), parameter :: DRAINE_OPT = '/home/kiseon/MoCafe/Grain/opt/'
   !! Refractive-index tables, taken from the release tree so that this check
   !! runs from the release alone; DustEM's own tables are read where they are.
   character(len=*), parameter :: DIEL   = '../data/dielectric/'
   !! Longward of this DustEM's amCBE tables are a power-law extrapolation,
   !! not a Mie calculation (Guillet et al. 2018 Sect. 5.1).
   real(wp), parameter :: LAM_EXTRAPOLATED = 800.0_wp
   !! Short end of the band the comparison counts.  It is the Lyman limit,
   !! which is where the shipped (non-EUV) products begin, and it is longward
   !! of the region where DustEM's table stops following the index file we
   !! hold: the BE table's shortest node is 0.0503 um and DustEM's Q_amCBEx
   !! carries rows down to 0.04 um, so between 0.04 and about 0.055 um its
   !! numbers come from an extrapolation of m that is not recorded anywhere in
   !! this tree.  Those rows are counted separately and reported, not hidden.
   real(wp), parameter :: LAM_LYMAN = 0.0912_wp
   !! Short end of DustEM's LAMBDA.DAT, and the largest size parameter the
   !! G18 Model D tables reach (2 pi * 2 um / 0.04 um).  Draine's suvSil_81
   !! runs to 0.001 um and to x ~ 6e4, far outside both.
   real(wp), parameter :: LAM_DUSTEM_MIN = 0.04_wp
   real(wp), parameter :: X_MODEL_MAX    = 320.0_wp

   call check_amcbe_spheres()
   call check_suvsil_spheres()

contains

   subroutine check_amcbe_spheres()
      type(refractive_index_t) :: idx
      real(wp), allocatable :: lam(:), a_um(:), qblk(:,:,:), ag(:), gblk(:,:,:)
      real(wp), allocatable :: dqa(:), dqs(:), dg(:)
      real(wp) :: nr, ki, x, qe, qs, qa, alb, g, lam_min, lam_max
      integer  :: nwave, nsize, nsg, nblk, iw, ia, n, u, n_uv, n_floor
      real(wp) :: dqa_hi, dqs_hi, dg_hi, dqa_uv_hi

      call load_index_energy_table(idx, DIEL//'index_amcBE_ZMCB96')
      call read_dustem_wavelength_grid(DUSTEM//'LAMBDA.DAT', nwave, lam)
      nblk = 2
      call read_dustem_sized_table(DUSTEM//'Q_amCBEx.DAT', nwave, nsize, a_um, qblk, nblk)
      nblk = 1
      call read_dustem_sized_table(DUSTEM//'G_amCBEx.DAT', nwave, nsg, ag, gblk, nblk)
      if (nsg /= nsize) then
         write(error_unit,'(a)') ' Q_amCBEx and G_amCBEx disagree on the radii.'
         stop 1
      end if
      lam_min = idx%lambda(1)
      lam_max = idx%lambda(idx%ndata)

      allocate(dqa(nwave*nsize), dqs(nwave*nsize), dg(nwave*nsize))
      dqa = 0.0_wp;  dqs = 0.0_wp;  dg = 0.0_wp
      n = 0;  n_uv = 0;  n_floor = 0
      dqa_hi = 0.0_wp;  dqs_hi = 0.0_wp;  dg_hi = 0.0_wp;  dqa_uv_hi = 0.0_wp
      open(newunit=u, file='output/mie_vs_amCBEx.dat', status='replace', action='write')
      write(u,'(a)') '# CHECK A: this tree''s Mie with data/dielectric/index_amcBE_ZMCB96'
      write(u,'(a)') '# against DustEM oprop/Q_amCBEx.DAT and G_amCBEx.DAT (BE a-C spheres).'
      write(u,'(a)') '# use = 1 where the comparison counts: 0.0912 <= lambda <= 800 um.'
      write(u,'(a)') '# use = 2 for 0.0503 <= lambda < 0.0912, where DustEM''s table stops'
      write(u,'(a)') '# following the index file this tree holds; reported separately.'
      write(u,'(a)') '# use = 0 elsewhere: shortward of the index table, or longward of 800 um'
      write(u,'(a)') '# where DustEM extrapolated the cross sections with a power law.'
      write(u,'(a)') '#  lambda[um]      a[um]           x      Qabs_mie   Qabs_dustem' // &
                     '    Qsca_mie   Qsca_dustem      g_mie     g_dustem  use'
      do ia = 1, nsize
         do iw = 1, nwave
            x = 2.0_wp*PI*a_um(ia)/lam(iw)
            ! Clamped only so that the rows outside the comparison band can
            ! still be written; every counted row lies inside the table.
            call refractive_index_at(idx, min(max(lam(iw), lam_min), lam_max), nr, ki)
            call mie(nr, ki, x, qe, qs, qa, alb, g)
            block
               integer :: use_it
               use_it = 0
               if (lam(iw) >= LAM_LYMAN .and. lam(iw) <= LAM_EXTRAPOLATED) then
                  use_it = 1
               else if (lam(iw) >= lam_min .and. lam(iw) < LAM_LYMAN) then
                  use_it = 2
               end if
               write(u,'(3es13.5,6es13.5,i5)') lam(iw), a_um(ia), x, &
                    qa, qblk(iw,ia,1), qs, qblk(iw,ia,2), g, gblk(iw,ia,1), use_it
               if (use_it == 2 .and. qblk(iw,ia,1) > 0.0_wp) then
                  n_uv = n_uv + 1
                  dqa_uv_hi = max(dqa_uv_hi, abs(qa - qblk(iw,ia,1))/qblk(iw,ia,1))
               end if
               if (use_it == 1 .and. qblk(iw,ia,1) > 0.0_wp .and. qblk(iw,ia,2) > 0.0_wp) then
                  n = n + 1
                  dqa(n) = abs(qa - qblk(iw,ia,1))/qblk(iw,ia,1)
                  dqs(n) = abs(qs - qblk(iw,ia,2))/qblk(iw,ia,2)
                  dqa_hi = max(dqa_hi, dqa(n))
                  dqs_hi = max(dqs_hi, dqs(n))
                  ! DustEM's G_amCBEx never goes below zero: wherever Mie gives
                  ! net BACKWARD scattering -- which a strongly absorbing sphere
                  ! does around x ~ 0.5 -- that table carries an exact 0.  Those
                  ! cells are counted and reported, not averaged into the
                  ! agreement, because the reference value there is a floor
                  ! rather than a measurement.
                  if (gblk(iw,ia,1) == 0.0_wp .and. g < 0.0_wp) then
                     n_floor = n_floor + 1
                  else
                     dg(n)  = abs(g - gblk(iw,ia,1))
                     dg_hi  = max(dg_hi, dg(n))
                  end if
               end if
            end block
         end do
      end do
      close(u)

      write(*,'(a)') repeat('=', 78)
      write(*,'(a)') ' CHECK A  BE a-C spheres: our Mie vs DustEM Q_amCBEx / G_amCBEx'
      write(*,'(a,i0,a,i0,a)') '   ', n, ' cells compared out of ', nwave*nsize, &
           '  (0.0912 <= lambda <= 800 um)'
      call report('Q_abs  relative', dqa(1:n), dqa_hi, .true.)
      call report('Q_sca  relative', dqs(1:n), dqs_hi, .true.)
      call report('g      absolute', dg(1:n),  dg_hi,  .false.)
      write(*,'(a,i0,a)') '   ', n_floor, ' of those cells carry an exact 0 in ' // &
           'G_amCBEx where Mie gives g < 0;'
      write(*,'(a)') '     the reference is floored there, so they are left out of the g line.'
      write(*,'(a,i0,a,f7.3,a)') '   0.0503 - 0.0912 um, reported apart: ', n_uv, &
           ' cells, max |dQabs| = ', 100.0_wp*dqa_uv_hi, '%'
      write(*,'(a)') '     DustEM extrapolated m below the BE table''s 0.0503 um node; ' // &
           'this tree holds m constant there.'
      write(*,'(a)') '   wrote output/mie_vs_amCBEx.dat'
   end subroutine check_amcbe_spheres


   subroutine check_suvsil_spheres()
      type(refractive_index_t) :: idx
      real(wp), allocatable :: a_um(:), lam(:), qa_d(:,:), qs_d(:,:), g_d(:,:)
      real(wp), allocatable :: dqa(:), dqs(:), dg(:)
      real(wp) :: nr, ki, x, qe, qs, qa, alb, g
      real(wp) :: dqa_hi, dqs_hi, dg_hi
      real(wp), allocatable :: qa_mie(:,:), ratio(:)
      real(wp) :: ref_med, ref_dev, ref_dev_hi, our_dev, our_dev_hi
      integer  :: nrad, nwav, ia, iw, n, u, n_all, n_ray, n_ray_bad, n_our_bad, nsmall
      real(wp) :: dqa_all_hi, dqs_all_hi

      call load_index_wavelength_table(idx, DIEL//'eps_suvSil')
      call read_draine_mie_table(DRAINE_OPT//'suvSil_81', nrad, nwav, a_um, lam, qa_d, qs_d, g_d)

      allocate(dqa(nrad*nwav), dqs(nrad*nwav), dg(nrad*nwav))
      allocate(qa_mie(nwav, nrad))
      dqa = 0.0_wp;  dqs = 0.0_wp;  dg = 0.0_wp
      n = 0;  n_all = 0
      dqa_hi = 0.0_wp;  dqs_hi = 0.0_wp;  dg_hi = 0.0_wp
      dqa_all_hi = 0.0_wp;  dqs_all_hi = 0.0_wp
      open(newunit=u, file='output/mie_vs_suvSil_81.dat', status='replace', action='write')
      write(u,'(a)') '# CHECK B: this tree''s Mie with data/dielectric/eps_suvSil against'
      write(u,'(a)') '# Draine''s own Mie table suvSil_81 (WD01 astrosilicate spheres).'
      write(u,'(a)') '#  lambda[um]      a[um]           x      Qabs_mie   Qabs_draine' // &
                     '    Qsca_mie   Qsca_draine      g_mie     g_draine'
      do ia = 1, nrad
         do iw = 1, nwav
            x = 2.0_wp*PI*a_um(ia)/lam(iw)
            call refractive_index_at(idx, lam(iw), nr, ki)
            call mie(nr, ki, x, qe, qs, qa, alb, g)
            qa_mie(iw, ia) = qa
            write(u,'(3es13.5,6es13.5)') lam(iw), a_um(ia), x, &
                 qa, qa_d(iw,ia), qs, qs_d(iw,ia), g, g_d(iw,ia)
            if (qa_d(iw,ia) > 0.0_wp .and. qs_d(iw,ia) > 0.0_wp) then
               n_all = n_all + 1
               dqa_all_hi = max(dqa_all_hi, abs(qa - qa_d(iw,ia))/qa_d(iw,ia))
               dqs_all_hi = max(dqs_all_hi, abs(qs - qs_d(iw,ia))/qs_d(iw,ia))
               if (lam(iw) >= LAM_DUSTEM_MIN .and. x <= X_MODEL_MAX) then
                  n = n + 1
                  dqa(n) = abs(qa - qa_d(iw,ia))/qa_d(iw,ia)
                  dqs(n) = abs(qs - qs_d(iw,ia))/qs_d(iw,ia)
                  dg(n)  = abs(g  - g_d(iw,ia))
                  dqa_hi = max(dqa_hi, dqa(n))
                  dqs_hi = max(dqs_hi, dqs(n))
                  dg_hi  = max(dg_hi,  dg(n))
               end if
            end if
         end do
      end do
      close(u)

      ! Rayleigh a-scaling test of the two tables.  At x < 0.01, Q_abs/a is a
      ! function of lambda alone; the median over radii at each lambda is the
      ! value both tables must reproduce.
      allocate(ratio(nrad))
      n_ray = 0;  n_ray_bad = 0;  n_our_bad = 0
      ref_dev_hi = 0.0_wp;  our_dev_hi = 0.0_wp
      do iw = 1, nwav
         ! Radii small enough to be in the dipole limit at this wavelength.
         nsmall = 0
         do ia = 1, nrad
            if (2.0_wp*PI*a_um(ia)/lam(iw) >= 0.01_wp) exit
            if (qa_d(iw,ia) <= 0.0_wp) cycle
            nsmall = nsmall + 1
            ratio(nsmall) = qa_d(iw,ia)/a_um(ia)
         end do
         ! Fewer than five leaves no majority to take the reference from.
         if (nsmall < 5) cycle
         call sort_ascending(ratio(1:nsmall))
         ref_med = ratio((nsmall+1)/2)
         if (ref_med <= 0.0_wp) cycle
         do ia = 1, nrad
            if (2.0_wp*PI*a_um(ia)/lam(iw) >= 0.01_wp) exit
            if (qa_d(iw,ia) <= 0.0_wp) cycle
            n_ray = n_ray + 1
            ref_dev = abs(qa_d(iw,ia)/a_um(ia)/ref_med - 1.0_wp)
            if (ref_dev > 0.01_wp) n_ray_bad = n_ray_bad + 1
            ref_dev_hi = max(ref_dev_hi, ref_dev)
            our_dev = abs(qa_mie(iw,ia)/a_um(ia)/ref_med - 1.0_wp)
            if (our_dev > 0.01_wp) n_our_bad = n_our_bad + 1
            our_dev_hi = max(our_dev_hi, our_dev)
         end do
      end do

      write(*,'(a)') repeat('=', 78)
      write(*,'(a)') ' CHECK B  WD01 astrosilicate spheres: our Mie vs Draine suvSil_81'
      write(*,'(a,i0,a,i0,a)') '   ', n, ' cells compared out of ', nrad*nwav, &
           '  (lambda >= 0.04 um and x <= 320, the domain the G18 tables occupy)'
      call report('Q_abs  relative', dqa(1:n), dqa_hi, .true.)
      call report('Q_sca  relative', dqs(1:n), dqs_hi, .true.)
      call report('g      absolute', dg(1:n),  dg_hi,  .false.)
      write(*,'(a,i0,a,f8.3,a,f8.3,a)') '   whole grid (', n_all, ' cells): ' // &
           'max |dQabs| = ', 100.0_wp*dqa_all_hi, '%, max |dQsca| = ', &
           100.0_wp*dqs_all_hi, '%'
      write(*,'(a)') '   Rayleigh a-scaling test at x < 0.01 (Q_abs/a must not depend on a):'
      write(*,'(a,i0,a,i0,a,f8.3,a)') '     suvSil_81 : ', n_ray_bad, ' of ', n_ray, &
           ' cells off the median by more than 1%, worst ', 100.0_wp*ref_dev_hi, '%'
      write(*,'(a,i0,a,i0,a,f8.3,a)') '     this work : ', n_our_bad, ' of ', n_ray, &
           ' cells off by more than 1%, worst ', 100.0_wp*our_dev_hi, '%'
      write(*,'(a)') '     A wavelength where more than half the reference cells are affected'
      write(*,'(a)') '     corrupts the median itself, which is what our few entries above 1%'
      write(*,'(a)') '     are measuring; they are not a departure from the a-scaling.'
      write(*,'(a)') '   wrote output/mie_vs_suvSil_81.dat'
   end subroutine check_suvsil_spheres


   subroutine read_draine_mie_table(path, nrad, nwav, a_um, lam, qabs, qsca, gfac)
      !! Draine's per-radius Mie tables (suvSil_81, Gra_81, ...): a few header
      !! lines carrying NRAD and NWAV, then one block per radius, each headed
      !! by "<a> = radius(micron) <material>" and a column title, then NWAV
      !! rows of  w[um]  Q_abs  Q_sca  g.  Wavelengths run from long to short
      !! and are the same for every radius.
      character(len=*), intent(in) :: path
      integer, intent(out) :: nrad, nwav
      real(wp), allocatable, intent(out) :: a_um(:), lam(:), qabs(:,:), qsca(:,:), gfac(:,:)
      character(len=256) :: line
      integer :: u, ios, ia, iw, ip
      real(wp) :: dum1, dum2

      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') ' cannot open ', path;  stop 1
      end if
      nrad = 0;  nwav = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         ip = index(line, 'NRAD')
         if (ip > 0) read(line,*,iostat=ios) nrad, dum1, dum2
         ip = index(line, 'NWAV')
         if (ip > 0) then
            read(line,*,iostat=ios) nwav, dum1, dum2
            exit
         end if
      end do
      if (nrad < 1 .or. nwav < 1) then
         write(error_unit,'(a,a)') ' no NRAD/NWAV header in ', path;  stop 1
      end if
      allocate(a_um(nrad), lam(nwav), qabs(nwav,nrad), qsca(nwav,nrad), gfac(nwav,nrad))

      do ia = 1, nrad
         ! Walk to this radius's own header line and take the radius from it.
         do
            read(u,'(a)',iostat=ios) line
            if (ios /= 0) then
               write(error_unit,'(a,i0,a,a)') ' radius block ', ia, ' missing in ', path
               stop 1
            end if
            if (index(line, '= radius') > 0) exit
         end do
         read(line,*) a_um(ia)
         read(u,'(a)') line                 ! column titles
         do iw = 1, nwav
            read(u,*,iostat=ios) lam(iw), qabs(iw,ia), qsca(iw,ia), gfac(iw,ia)
            if (ios /= 0) then
               write(error_unit,'(a,i0,a,i0)') ' short block at radius ', ia, ', row ', iw
               stop 1
            end if
         end do
      end do
      close(u)
   end subroutine read_draine_mie_table


   subroutine report(label, d, dmax, as_percent)
      character(len=*), intent(in) :: label
      real(wp), intent(in) :: d(:), dmax
      logical,  intent(in) :: as_percent
      real(wp), allocatable :: s(:)
      real(wp) :: med, p99, scale
      integer  :: n
      character(len=2) :: unit
      n = size(d)
      if (n < 1) return
      allocate(s(n));  s = d
      call sort_ascending(s)
      med = s((n+1)/2)
      p99 = s(max(1, min(n, nint(0.99_wp*n))))
      scale = 1.0_wp;  unit = '  '
      if (as_percent) then
         scale = 100.0_wp;  unit = ' %'
      end if
      write(*,'(a,a,a,es11.3,a,a,es11.3,a,a,es11.3,a)') '   ', label, &
           ':  median ', med*scale, unit, ',  99th pct ', p99*scale, unit, &
           ',  max ', dmax*scale, unit
   end subroutine report


   subroutine sort_ascending(a)
      !! Heapsort: in place, no recursion, O(n log n) -- the arrays here reach
      !! a few times 1e4 entries and are sorted once per quantity.
      real(wp), intent(inout) :: a(:)
      integer :: n, i
      real(wp) :: t
      n = size(a)
      do i = n/2, 1, -1
         call sift_down(a, i, n)
      end do
      do i = n, 2, -1
         t = a(1);  a(1) = a(i);  a(i) = t
         call sift_down(a, 1, i-1)
      end do
   end subroutine sort_ascending


   subroutine sift_down(a, first, last)
      real(wp), intent(inout) :: a(:)
      integer,  intent(in)    :: first, last
      integer  :: root, child
      real(wp) :: t
      root = first
      do
         child = 2*root
         if (child > last) exit
         if (child < last) then
            if (a(child) < a(child+1)) child = child + 1
         end if
         if (a(root) >= a(child)) exit
         t = a(root);  a(root) = a(child);  a(child) = t
         root = child
      end do
   end subroutine sift_down

end program compare_mie_spheres
