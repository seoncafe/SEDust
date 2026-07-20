program compare_go_jori
   ! Geometric-optics check of the orientation-resolved cross sections.
   !
   ! Over a strided subset of the geometric-optics regime (x > 50, with
   ! x = 2 pi a_eff / lambda) of the DH21 astrodust grid, this computes the
   ! oriented efficiencies Q_ext(jori), Q_abs(jori), Q_sca(jori) from
   ! geometric_optics_limit (projected-area extinction paradox for Q_ext, a
   ! surface Fresnel integral for the opaque-grain Q_abs) and compares the
   ! derived combinations, in the convention of sed/src/q_table_jori.f90:
   !
   !   qpol = 0.5 * (Q(jori=3) - Q(jori=2))     polarized part
   !   qran = (Q(1) + Q(2) + Q(3)) / 3          orientation mean
   !
   ! against the Q_ext and Q_abs blocks of the Hensley & Draine (2021/2023)
   ! orientation-resolved table.
   !
   ! Anchors reported:
   !   phys  physical bounds at every node: Q_abs(jori) in [0, Q_ext(jori)],
   !         Q_sca(jori) >= 0, albedo in [0,1].  Violations counted.
   !   E     extinction continuity across x = 50: the T-matrix oriented
   !         Q_ext(jori) at x just below 50 must connect to the projected-area
   !         Q_ext(jori) at x just above 50.
   !   A     absorption continuity across x = 50 (the decisive one): the
   !         T-matrix oriented Q_abs(jori) and qpol_abs at x just below 50 must
   !         connect to the Fresnel-GO Q_abs(jori) and qpol_abs at x just above.
   !   conv  Fresnel surface-integral convergence: recompute Q_abs(jori) with
   !         the quadrature doubled and report the relative change.
   !
   ! The refractive index m(lambda), the spheroid shape (oblate b/a = 1.4), and
   ! the x-window match run_tmatrix.f90.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use constants,        only: wp
   use read_index,       only: load_index, interp_m
   use asymptotic_optics, only: geometric_optics_limit, &
                                projected_area_extinction, fresnel_opaque_absorption
   use tmatrix_oriented, only: tmatrix_oriented_cross
   implicit none

   character(len=*), parameter :: f_aeff  = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: f_wave  = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'
   character(len=*), parameter :: f_qjori = &
      '../data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: f_scratch = 'q_jori_go_cmp_scratch.dat'

   integer,  parameter :: NA = 169, NW = 1129, NORI = 3, NHEAD = 12
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2, NP_OBL = -1
   real(wp), parameter :: X_LARGE = 50.0_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   ! Main comparison stride over the x > 50 region.  GO is cheap (no T-matrix
   ! solve), so a moderate stride still spans the region well.
   integer,  parameter :: JW_STEP = 4, JA_STEP = 2

   real(wp) :: a_eff(NA), lambda(NW)
   real(wp) :: nr_cache(NW), ki_cache(NW)
   real(wp), allocatable :: qext_hd(:,:,:), qabs_hd(:,:,:)   ! (jw, ja, jori)
   real(wp), allocatable :: dev_ranE(:), dev_polE(:), dev_ranA(:), dev_polA(:)
   real(wp) :: qe(NORI), qa(NORI), qs(NORI)
   real(wp) :: qe_avg, qs_avg, walb_avg, g_avg
   real(wp) :: x, nr, ki, w
   real(wp) :: qpolE_o, qranE_o, qpolA_o, qranA_o
   real(wp) :: qpolE_h, qranE_h, qpolA_h, qranA_h, floor_pol
   integer  :: jw, ja, jo, n_win, n_ranE, n_polE, n_ranA, n_polA, n_viol
   real(wp) :: t0, t1

   call cpu_time(t0)

   call read_one_col(f_aeff, NA, a_eff)
   call read_one_col(f_wave, NW, lambda)
   call load_index(f_index)
   do jw = 1, NW
      call interp_m(lambda(jw), nr_cache(jw), ki_cache(jw))
   end do

   allocate(qext_hd(NW, NA, NORI), qabs_hd(NW, NA, NORI))
   call read_q_blocks(f_qjori, f_scratch, qext_hd, qabs_hd)

   ! ---- Anchor E and A: continuity across x = 50 -----------------------
   call anchor_continuity()

   ! ---- Anchor conv: Fresnel surface-integral convergence --------------
   call anchor_quadrature_convergence()

   ! ---- main comparison sweep over x > 50 ------------------------------
   allocate(dev_ranE(NW*NA), dev_polE(NW*NA), dev_ranA(NW*NA), dev_polA(NW*NA))
   n_win = 0;  n_ranE = 0;  n_polE = 0;  n_ranA = 0;  n_polA = 0;  n_viol = 0

   do jw = 1, NW, JW_STEP
      nr = nr_cache(jw);  ki = ki_cache(jw)
      do ja = 1, NA, JA_STEP
         x = 2.0_wp * PI * a_eff(ja) / lambda(jw)
         if (x <= X_LARGE) cycle
         n_win = n_win + 1

         call geometric_optics_limit(a_eff(ja), lambda(jw), nr, ki, EPS_BA, &
                                     qe_avg, qs_avg, walb_avg, g_avg, &
                                     qext_ori=qe, qabs_ori=qa, qsca_ori=qs)

         ! Anchor phys: physical bounds per orientation.
         do jo = 1, NORI
            w = 0.0_wp
            if (qe(jo) > 0.0_wp) w = qs(jo) / qe(jo)
            if (qa(jo) < 0.0_wp .or. qa(jo) > qe(jo) + 1.0e-9_wp .or. &
                qs(jo) < 0.0_wp .or. w < -1.0e-9_wp .or. w > 1.0_wp + 1.0e-9_wp) then
               n_viol = n_viol + 1
               write(error_unit,'(a,i0,a,3(1x,es12.4))') &
                  ' PHYS VIOLATION jori=', jo, '  Qext,Qabs,Qsca =', &
                  qe(jo), qa(jo), qs(jo)
               write(error_unit,'(a,es12.4,a,es12.4,a,es12.4)') &
                  '   at lambda=', lambda(jw), '  a_eff=', a_eff(ja), '  x=', x
            end if
         end do

         qranE_o = (qe(1) + qe(2) + qe(3)) / 3.0_wp
         qpolE_o = 0.5_wp * (qe(3) - qe(2))
         qranA_o = (qa(1) + qa(2) + qa(3)) / 3.0_wp
         qpolA_o = 0.5_wp * (qa(3) - qa(2))

         qranE_h = (qext_hd(jw,ja,1) + qext_hd(jw,ja,2) + qext_hd(jw,ja,3)) / 3.0_wp
         qpolE_h = 0.5_wp * (qext_hd(jw,ja,3) - qext_hd(jw,ja,2))
         qranA_h = (qabs_hd(jw,ja,1) + qabs_hd(jw,ja,2) + qabs_hd(jw,ja,3)) / 3.0_wp
         qpolA_h = 0.5_wp * (qabs_hd(jw,ja,3) - qabs_hd(jw,ja,2))

         if (abs(qranE_h) > 0.0_wp) then
            n_ranE = n_ranE + 1;  dev_ranE(n_ranE) = abs(qranE_o/qranE_h - 1.0_wp)
         end if
         if (abs(qranA_h) > 0.0_wp) then
            n_ranA = n_ranA + 1;  dev_ranA(n_ranA) = abs(qranA_o/qranA_h - 1.0_wp)
         end if
         ! Polarized parts can pass through zero: ratio only where the HD23
         ! signal is meaningful relative to the local mean.
         floor_pol = 1.0e-4_wp * abs(qranE_h)
         if (abs(qpolE_h) > floor_pol) then
            n_polE = n_polE + 1;  dev_polE(n_polE) = abs(qpolE_o/qpolE_h - 1.0_wp)
         end if
         floor_pol = 1.0e-4_wp * abs(qranA_h)
         if (abs(qpolA_h) > floor_pol) then
            n_polA = n_polA + 1;  dev_polA(n_polA) = abs(qpolA_o/qpolA_h - 1.0_wp)
         end if
      end do
   end do

   call cpu_time(t1)

   write(*,'(a)') '======================================================================'
   write(*,'(a)') ' Geometric-optics oriented Q vs HD23 (x > 50)'
   write(*,'(a)') '======================================================================'
   write(*,'(a,i0,a,i0)') ' lambda stride : every ', JW_STEP, ' of ', NW
   write(*,'(a,i0,a,i0)') ' a_eff  stride : every ', JA_STEP, ' of ', NA
   write(*,'(a,i0)')      ' strided nodes with x > 50 : ', n_win
   write(*,'(a,i0)')      ' phys-bound violations     : ', n_viol
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' qran_ext = (Q1+Q2+Q3)/3   [extinction paradox + projected-area split]'
   write(*,'(a,i0)')     '   nodes                  : ', n_ranE
   write(*,'(a,es12.4)') '   median |ours/HD23 - 1| : ', median(dev_ranE, n_ranE)
   write(*,'(a,es12.4)') '   max    |ours/HD23 - 1| : ', maxval_n(dev_ranE, n_ranE)
   write(*,'(a)') ' qpol_ext = 0.5*(Q3-Q2)    [ours = 0 at this order]'
   write(*,'(a,i0,a)')   '   nodes (|qpolE_HD|>1e-4 qranE_HD) : ', n_polE, &
        '   -- ours = 0, so |ours/HD - 1| = 1 by construction'
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a)') ' qran_abs = (Q1+Q2+Q3)/3   [Fresnel opaque absorption]'
   write(*,'(a,i0)')     '   nodes                  : ', n_ranA
   write(*,'(a,es12.4)') '   median |ours/HD23 - 1| : ', median(dev_ranA, n_ranA)
   write(*,'(a,es12.4)') '   max    |ours/HD23 - 1| : ', maxval_n(dev_ranA, n_ranA)
   write(*,'(a)') ' qpol_abs = 0.5*(Q3-Q2)    [absorption dichroism]'
   write(*,'(a,i0,a)')   '   nodes (|qpolA_HD|>1e-4 qranA_HD) : ', n_polA
   write(*,'(a,es12.4)') '   median |ours/HD23 - 1| : ', median(dev_polA, n_polA)
   write(*,'(a,es12.4)') '   max    |ours/HD23 - 1| : ', maxval_n(dev_polA, n_polA)
   write(*,'(a)') '----------------------------------------------------------------------'
   write(*,'(a,f8.1,a)') ' wall (cpu) time : ', t1 - t0, ' s'
   write(*,'(a)') '======================================================================'

contains

   subroutine anchor_continuity()
      ! For a handful of large radii, find the grid node whose x is just below
      ! 50 (T-matrix oriented cross section) and the node just above 50
      ! (geometric-optics oriented cross section), and print Q_ext(jori),
      ! Q_abs(jori), qpol_ext, qpol_abs on both sides so the connection across
      ! x = 50 can be judged.  If the T-matrix solve nearest below 50 does not
      ! converge, step to smaller x until it does (the x used is printed).
      ! Each T-matrix solve near x = 50 costs ~50 s, so only a few of the
      ! largest radii are used and the below-50 search is capped.
      integer, parameter :: NSHOW = 3, MAX_TRY = 4
      integer  :: pick(NSHOW), kk, ja, jw, jw_lo, jw_hi, kerr, ntry
      real(wp) :: xlo, xhi, mm_dummy
      real(wp) :: qe_t(NORI), qs_t(NORI), qa_t(NORI)
      real(wp) :: qe_g(NORI), qs_g(NORI), qa_g(NORI)
      real(wp) :: qe_avg, qs_avg, wa_avg, g_avg, xt
      complex(wp) :: mm

      ! Largest five radii on the grid (they reach x > 50 at the UV end).
      do kk = 1, NSHOW
         pick(kk) = NA - NSHOW + kk
      end do

      write(*,'(a)') '======================================================================'
      write(*,'(a)') ' Anchor E & A: continuity across x = 50'
      write(*,'(a)') '   T-matrix oriented (x just below 50) | GO oriented (x just above 50)'
      write(*,'(a)') '======================================================================'

      do kk = 1, NSHOW
         ja = pick(kk)
         ! GO side: smallest x > 50 -> largest lambda giving x > 50.  Lambda
         ! ascends, x = 2 pi a / lambda descends, so scan lambda upward and
         ! take the last node still above 50.
         jw_hi = 0
         do jw = 1, NW
            xhi = 2.0_wp*PI*a_eff(ja)/lambda(jw)
            if (xhi > X_LARGE) jw_hi = jw
         end do
         if (jw_hi == 0) cycle
         xhi = 2.0_wp*PI*a_eff(ja)/lambda(jw_hi)
         call geometric_optics_limit(a_eff(ja), lambda(jw_hi), &
                 nr_cache(jw_hi), ki_cache(jw_hi), EPS_BA, &
                 qe_avg, qs_avg, wa_avg, g_avg, &
                 qext_ori=qe_g, qabs_ori=qa_g, qsca_ori=qs_g)

         ! T-matrix side: smallest lambda with x < 50 that converges (capped).
         jw_lo = jw_hi + 1
         kerr = 1;  ntry = 0
         do jw = jw_lo, NW
            ntry = ntry + 1
            if (ntry > MAX_TRY) exit
            xt = 2.0_wp*PI*a_eff(ja)/lambda(jw)
            mm = cmplx(nr_cache(jw), ki_cache(jw), kind=wp)
            call tmatrix_oriented_cross(a_eff(ja), lambda(jw), mm, EPS_BA, &
                                        NP_OBL, DDELT, NDGS, qe_t, qs_t, qa_t, kerr)
            if (kerr == 0) then
               jw_lo = jw;  xlo = xt;  exit
            end if
         end do

         write(*,'(a)') '----------------------------------------------------------------------'
         write(*,'(a,es11.4,a)') ' a_eff = ', a_eff(ja), ' um'
         if (kerr == 0) then
            write(*,'(a,f8.3,a,es12.4,es12.4,es12.4)') &
               '   Tmat  x=', xlo, '  Qext(1,2,3)=', qe_t(1), qe_t(2), qe_t(3)
            write(*,'(a,es12.4,es12.4,es12.4)') &
               '                    Qabs(1,2,3)=', qa_t(1), qa_t(2), qa_t(3)
            write(*,'(a,es12.4,a,es12.4)') &
               '                    qpolE=', 0.5_wp*(qe_t(3)-qe_t(2)), &
               '  qpolA=', 0.5_wp*(qa_t(3)-qa_t(2))
         else
            write(*,'(a)') '   Tmat  (no converged node with x < 50)'
         end if
         write(*,'(a,f8.3,a,es12.4,es12.4,es12.4)') &
            '   GO    x=', xhi, '  Qext(1,2,3)=', qe_g(1), qe_g(2), qe_g(3)
         write(*,'(a,es12.4,es12.4,es12.4)') &
            '                    Qabs(1,2,3)=', qa_g(1), qa_g(2), qa_g(3)
         write(*,'(a,es12.4,a,es12.4)') &
            '                    qpolE=', 0.5_wp*(qe_g(3)-qe_g(2)), &
            '  qpolA=', 0.5_wp*(qa_g(3)-qa_g(2))
      end do
      write(*,'(a)') '======================================================================'
      mm_dummy = 0.0_wp
   end subroutine anchor_continuity

   subroutine anchor_quadrature_convergence()
      ! Recompute the Fresnel surface integral with the quadrature doubled at
      ! representative (m, x) points spanning the UV index range, and report
      ! the relative change in Q_abs(jori=2) and Q_abs(jori=3).
      integer, parameter :: NPT = 3
      integer  :: k, jw
      real(wp) :: e_targets(NPT), best
      integer  :: jw_pick(NPT)
      real(wp) :: qa2_1, qa2_2, qa3_1, qa3_2, ap
      complex(wp) :: mm
      real(wp), parameter :: KHAT_PERP(3) = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
      real(wp), parameter :: EHAT_AXIS(3) = (/ 0.0_wp, 0.0_wp, 1.0_wp /)
      real(wp), parameter :: EHAT_EQ(3)   = (/ 0.0_wp, 1.0_wp, 0.0_wp /)

      ! Target photon energies [eV] in the UV band that x > 50 samples.
      e_targets = (/ 5.0_wp, 9.0_wp, 12.0_wp /)
      ! DH21_wave is in um; nearest-energy node via lambda_um = 1.23984/E_eV.
      do k = 1, NPT
         best = huge(1.0_wp);  jw_pick(k) = 1
         do jw = 1, NW
            if (abs(1.23984_wp/lambda(jw) - e_targets(k)) < best) then
               best = abs(1.23984_wp/lambda(jw) - e_targets(k));  jw_pick(k) = jw
            end if
         end do
      end do

      write(*,'(a)') '======================================================================'
      write(*,'(a)') ' Anchor conv: Fresnel surface-integral convergence (double quadrature)'
      write(*,'(a)') '======================================================================'
      write(*,'(a)') '     E[eV]     n_r      k_i    |Qabs2(2x)/-1|   |Qabs3(2x)/-1|'
      do k = 1, NPT
         jw = jw_pick(k)
         mm = cmplx(nr_cache(jw), ki_cache(jw), kind=wp)
         call fresnel_opaque_absorption(EPS_BA, mm, KHAT_PERP, EHAT_AXIS, 64, qa2_1, ap)
         call fresnel_opaque_absorption(EPS_BA, mm, KHAT_PERP, EHAT_AXIS,128, qa2_2)
         call fresnel_opaque_absorption(EPS_BA, mm, KHAT_PERP, EHAT_EQ,   64, qa3_1)
         call fresnel_opaque_absorption(EPS_BA, mm, KHAT_PERP, EHAT_EQ,  128, qa3_2)
         write(*,'(f9.2,f9.4,f9.4,3x,es13.4,3x,es13.4)') &
            1.23984_wp/lambda(jw), nr_cache(jw), ki_cache(jw), &
            abs(qa2_2/qa2_1 - 1.0_wp), abs(qa3_2/qa3_1 - 1.0_wp)
      end do
      ! Also confirm the numeric projected area matches pi a_s c (jori=2/3).
      write(*,'(a,es13.5,a,es13.5)') '   numeric A_proj(jori=2/3) = ', ap, &
           '   analytic pi a_s c = ', PI*EPS_BA**(-1.0_wp/3.0_wp)
      write(*,'(a)') '======================================================================'
   end subroutine anchor_quadrature_convergence

   subroutine read_one_col(filename, n, arr)
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: n
      real(wp),         intent(out) :: arr(n)
      integer :: u, ios
      character(len=512) :: header
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') ' ERROR: cannot open ', trim(filename)
         stop 1
      end if
      read(u,'(a)') header
      read(u,'(a)') header
      read(u,*) arr(1:n)
      close(u)
   end subroutine read_one_col

   subroutine read_q_blocks(gz_file, scratch, qe_tab, qa_tab)
      ! Decompress the gzip'd HD23 table once and read its Q_ext (iq=1) and
      ! Q_abs (iq=2) blocks; skip Q_sca (iq=3).  On-disk stream: 12 header
      ! lines, then for each quantity, for each orientation jori = 1..3, for
      ! each wavelength (NW records), one record of the NA radii.
      character(len=*), intent(in)  :: gz_file, scratch
      real(wp),         intent(out) :: qe_tab(NW,NA,NORI), qa_tab(NW,NA,NORI)
      integer :: u, ios, estat, cstat, iq, jori, jwl, i
      real(wp) :: row(NA)
      character(len=512) :: line
      call execute_command_line('gzip -dc "'//trim(gz_file)//'" > "'// &
                                trim(scratch)//'"', exitstat=estat, cmdstat=cstat)
      if (cstat /= 0 .or. estat /= 0) then
         write(error_unit,'(a,a)') ' ERROR: gzip -dc failed on ', trim(gz_file)
         call delete_scratch(scratch);  stop 1
      end if
      open(newunit=u, file=trim(scratch), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') ' ERROR: cannot open ', trim(scratch)
         call delete_scratch(scratch);  stop 1
      end if
      do i = 1, NHEAD
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then
            write(error_unit,'(a)') ' ERROR: unexpected EOF in header'
            close(u);  call delete_scratch(scratch);  stop 1
         end if
      end do
      do iq = 1, 3
         do jori = 1, NORI
            do jwl = 1, NW
               read(u,*,iostat=ios) row(1:NA)
               if (ios /= 0) then
                  write(error_unit,'(a)') ' ERROR: read error in Q block'
                  close(u);  call delete_scratch(scratch);  stop 1
               end if
               if (iq == 1) qe_tab(jwl,:,jori) = row(1:NA)
               if (iq == 2) qa_tab(jwl,:,jori) = row(1:NA)
            end do
         end do
      end do
      close(u)
      call delete_scratch(scratch)
   end subroutine read_q_blocks

   subroutine delete_scratch(path)
      character(len=*), intent(in) :: path
      integer :: u, ios
      open(newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close(u, status='delete')
   end subroutine delete_scratch

   real(wp) function median(v, n)
      real(wp), intent(in) :: v(:)
      integer,  intent(in) :: n
      real(wp), allocatable :: s(:)
      integer :: i, j
      real(wp) :: t
      if (n <= 0) then
         median = 0.0_wp;  return
      end if
      allocate(s(n));  s = v(1:n)
      do i = 2, n
         t = s(i);  j = i - 1
         do while (j >= 1)
            if (s(j) <= t) exit
            s(j+1) = s(j);  j = j - 1
         end do
         s(j+1) = t
      end do
      if (mod(n,2) == 1) then
         median = s((n+1)/2)
      else
         median = 0.5_wp * (s(n/2) + s(n/2 + 1))
      end if
      deallocate(s)
   end function median

   real(wp) function maxval_n(v, n)
      real(wp), intent(in) :: v(:)
      integer,  intent(in) :: n
      if (n <= 0) then
         maxval_n = 0.0_wp
      else
         maxval_n = maxval(v(1:n))
      end if
   end function maxval_n

end program compare_go_jori
