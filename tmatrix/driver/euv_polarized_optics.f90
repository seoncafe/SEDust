program euv_polarized_optics
   ! Orientation-resolved optics of the DH21 astrodust spheroid (b/a = 1.4) in
   ! the extreme ultraviolet, i.e. SHORTWARD of the 0.0912 um (13.6 eV) short-
   ! wavelength end of the tabulated grid, computed from first principles at
   ! arbitrary (lambda, a_eff) nodes.
   !
   ! WHAT IT MEASURES.  For each node it evaluates, independently, all three
   ! descriptions of the same spheroid that this tree owns:
   !
   !   R  rayleigh_limit          dipole polarizability, exact as x -> 0
   !   T  T-matrix                Mishchenko, forward-amplitude optical theorem
   !                              (extinction and its real-part twin) plus the
   !                              fixed-orientation scattering-sphere integral
   !   G  geometric_optics_limit  extinction paradox + Fresnel surface integral
   !
   ! and prints them side by side, so the size-parameter range over which each
   ! is trustworthy -- and whether they agree where they overlap -- is read off
   ! the output rather than assumed.  The T-matrix status IERR is reported per
   ! node and NEVER hidden: a non-converged node is written with its Q left at
   ! zero and ierr /= 0, so no unconverged number can be mistaken for a result.
   !
   ! Polarized combinations formed from the orientation-resolved Q, in the
   ! convention of sed/src/q_table_jori.f90 (a = spheroid symmetry axis;
   ! jori=1: k || a, jori=2: k perp a and E || a, jori=3: k perp a and E perp a):
   !
   !   Q_pol,ext = 0.5*(Q_ext(3) - Q_ext(2))      dichroic extinction
   !   Q_bir,ext = 0.5*(Q_re (3) - Q_re (2))      birefringence
   !   Q_pol,abs = 0.5*(Q_abs(3) - Q_abs(2))      dichroic absorption
   !   Q_ran     = (Q(1) + Q(2) + Q(3))/3         random-orientation trace mean
   !
   ! The volume-equivalent Mie SPHERE is evaluated at the same node as well,
   ! since that is the optics the SED model currently uses in this band.
   !
   ! Usage (run from tmatrix/):
   !   ./euv_polarized_optics.x NODEFILE OUTFILE [sca|nosca] [XTM_HI]
   !
   !   NODEFILE  text, one "lambda[um]  a_eff[um]" pair per line; # comments ok
   !   OUTFILE   one output record per node (columns listed in its header)
   !   sca       also run the scattering-sphere integral, giving Q_sca and
   !             hence Q_abs = Q_ext - Q_sca per orientation (default);
   !             nosca skips it and leaves Q_sca / Q_abs at zero, which costs
   !             two AMPL calls instead of ~2*NMAX^2 and is enough when only
   !             the extinction and birefringence channels are wanted
   !   XTM_HI    largest x at which the T-matrix is attempted (default 200);
   !             the solver's own storage limit NPN1 = 100 caps the multipole
   !             order, so failures above x ~ 80 are expected and are exactly
   !             what this driver is meant to locate
   !   DDELT     T-matrix convergence tolerance (default 1e-3, the production
   !             value); tighten it to test whether a converged node is real
   !   NDGS      Gauss-quadrature multiplier (default 2, the production value)

   use, intrinsic :: iso_fortran_env, only: int64, error_unit
   use tmatrix_api,        only: wp, tmatrix_options_t, tmatrix_result_t, &
                                 tmatrix_scatmat_t, tmatrix_workspace_t, &
                                 tmatrix_workspace_init, tmatrix_workspace_finalize, &
                                 tmatrix_eval_scattering_matrix, ampl
   use tmatrix_status,     only: TMATRIX_SUCCESS, TMATRIX_ERR_NUMERICAL_NONFINITE
   use read_index,         only: load_index, interp_m
   use asymptotic_optics,  only: rayleigh_limit, geometric_optics_limit
   use tmatrix_oriented,   only: tmatrix_oriented_cross
   use mie_mod,            only: mie
   implicit none

   character(len=*), parameter :: F_INDEX = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'

   real(wp), parameter :: EPS_BA = 1.4_wp
   real(wp), parameter :: DDELT_DEF = 1.0e-3_wp
   integer,  parameter :: NDGS_DEF  = 2
   integer,  parameter :: NP_OBL = -1
   real(wp), parameter :: PI     = acos(-1.0_wp)
   ! Smallest x at which a T-matrix solve is attempted.  Below it the dipole
   ! limit is exact and the solver is both unnecessary and ill-conditioned.
   real(wp), parameter :: XTM_LO = 0.05_wp

   character(len=512) :: nodefile, outfile, arg
   real(wp), allocatable :: lam_in(:), a_in(:)
   integer  :: nnode, i, uo, ios
   logical  :: do_sca
   real(wp) :: xtm_hi, ddelt
   integer  :: ndgs

   real(wp) :: lam, a, x, nr, ki, mm1
   complex(wp) :: m
   real(wp) :: qe_r(3), qa_r(3), qs_r(3), qre_r(3)
   real(wp) :: qe_g(3), qa_g(3), qs_g(3), qre_g(3)
   real(wp) :: qe_t(3), qa_t(3), qs_t(3), qre_t(3)
   real(wp) :: qext_av, qsca_av, walb, asymm
   real(wp) :: qext_mie, qsca_mie, qabs_mie, alb_mie, g_mie
   integer  :: ierr
   integer(int64) :: c0, c1, crate
   real(wp) :: dt

   ! One workspace for the whole node list: each solve is serial and the two
   ! forward AMPL evaluations read the T-matrix it left behind.
   type(tmatrix_workspace_t) :: work
   integer  :: tm_status
   character(len=160) :: tm_message

   ! ---- command line ------------------------------------------------------
   if (command_argument_count() < 2) then
      write(*,'(a)') ' usage: euv_polarized_optics.x NODEFILE OUTFILE [sca|nosca] [XTM_HI]'
      stop 1
   end if
   call get_command_argument(1, nodefile)
   call get_command_argument(2, outfile)
   do_sca = .true.
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg)
      if (trim(arg) == 'nosca') do_sca = .false.
   end if
   xtm_hi = 200.0_wp
   if (command_argument_count() >= 4) then
      call get_command_argument(4, arg)
      read(arg,*,iostat=ios) xtm_hi
      if (ios /= 0) then
         write(error_unit,'(a)') ' bad XTM_HI'
         stop 1
      end if
   end if

   ddelt = DDELT_DEF
   if (command_argument_count() >= 5) then
      call get_command_argument(5, arg)
      read(arg,*,iostat=ios) ddelt
      if (ios /= 0) then
         write(error_unit,'(a)') ' bad DDELT'
         stop 1
      end if
   end if
   ndgs = NDGS_DEF
   if (command_argument_count() >= 6) then
      call get_command_argument(6, arg)
      read(arg,*,iostat=ios) ndgs
      if (ios /= 0) then
         write(error_unit,'(a)') ' bad NDGS'
         stop 1
      end if
   end if

   call tmatrix_workspace_init(work, tm_status, tm_message)
   if (tm_status /= TMATRIX_SUCCESS) then
      write(error_unit,'(a,a)') ' ERROR: libtmatrix: ', trim(tm_message)
      stop 2
   end if

   call read_nodes(trim(nodefile), lam_in, a_in, nnode)
   call load_index(F_INDEX)
   call system_clock(count_rate=crate)

   open(newunit=uo, file=trim(outfile), status='replace', action='write')
   call write_header(uo, do_sca, xtm_hi, ddelt, ndgs)

   do i = 1, nnode
      lam = lam_in(i)
      a   = a_in(i)
      x   = 2.0_wp*PI*a/lam
      call interp_m(lam, nr, ki)
      m   = cmplx(nr, ki, kind=wp)
      mm1 = abs(m - cmplx(1.0_wp, 0.0_wp, kind=wp))

      ! Rayleigh dipole limit -- always evaluated, exact only for x << 1.
      call rayleigh_limit(a, lam, nr, ki, EPS_BA, qext_av, qsca_av, walb, asymm, &
                          qext_ori=qe_r, qabs_ori=qa_r, qsca_ori=qs_r, qre_ori=qre_r)

      ! Geometric-optics limit -- always evaluated, exact only for x >> 1.
      call geometric_optics_limit(a, lam, nr, ki, EPS_BA, qext_av, qsca_av, walb, asymm, &
                          qext_ori=qe_g, qabs_ori=qa_g, qsca_ori=qs_g, qre_ori=qre_g)

      ! T-matrix -- attempted inside the x window; ierr reported verbatim.
      qe_t = 0.0_wp;  qa_t = 0.0_wp;  qs_t = 0.0_wp;  qre_t = 0.0_wp
      ierr = -1
      dt   = 0.0_wp
      if (x >= XTM_LO .and. x <= xtm_hi) then
         call system_clock(c0)
         if (do_sca) then
            call tmatrix_oriented_cross(work, a, lam, m, EPS_BA, NP_OBL, ddelt, ndgs, &
                                        qe_t, qs_t, qa_t, ierr, qre_ori=qre_t)
         else
            call forward_amplitude_extinction(work, a, lam, m, EPS_BA, NP_OBL, ddelt, ndgs, &
                                              qe_t, qre_t, ierr)
         end if
         call system_clock(c1)
         dt = real(c1 - c0, kind=wp) / real(crate, kind=wp)
         if (ierr /= 0) then
            qe_t = 0.0_wp;  qa_t = 0.0_wp;  qs_t = 0.0_wp;  qre_t = 0.0_wp
         end if
      end if

      ! Volume-equivalent Mie sphere: the optics the SED model uses here.
      call mie(nr, ki, x, qext_mie, qsca_mie, qabs_mie, alb_mie, g_mie)

      write(uo,'(2es14.6, es12.4, 2f10.5, es12.4, i5, es11.3, 3(3es14.6), &
                &3(3es14.6), 3(3es14.6), 5es14.6)') &
         lam, a, x, nr, ki, mm1, ierr, dt, &
         qe_r, qre_r, qa_r, &
         qe_t, qre_t, qa_t, &
         qe_g, qre_g, qa_g, &
         qext_mie, qabs_mie, qs_r(1), qs_t(1), qs_g(1)
      flush(uo)

      write(*,'(a,i0,a,i0,a,es11.4,a,es11.4,a,es10.3,a,i3,a,f8.2,a)') &
         ' node ', i, '/', nnode, '  lam=', lam, ' a=', a, ' x=', x, &
         ' ierr=', ierr, '  (', dt, ' s)'
   end do
   close(uo)
   call tmatrix_workspace_finalize(work)
   write(*,'(a,a)') ' wrote ', trim(outfile)

contains

   subroutine write_header(u, sca, xhi, dd, ng)
      integer,  intent(in) :: u, ng
      logical,  intent(in) :: sca
      real(wp), intent(in) :: xhi, dd
      write(u,'(a)') '# DH21 astrodust spheroid b/a = 1.4, orientation-resolved EUV optics'
      write(u,'(a)') '# jori=1: k || a ; jori=2: k perp a, E || a ; jori=3: k perp a, E perp a'
      write(u,'(a,l1,a,es10.3,a,es10.3,a,i0)') '# scattering integral = ', sca, &
         '   XTM_HI = ', xhi, '   DDELT = ', dd, '   NDGS = ', ng
      write(u,'(a)') '# ierr: -1 T-matrix not attempted (x outside window), 0 converged,'
      write(u,'(a)') '#       >0 T-matrix status (1 NMAX>NPN1, 2 NGAUSS>NPNG1, 3 no'
      write(u,'(a)') '#       convergence at NPN1, 4 NMAX1>NPN4, 5 Gauss refinement).'
      write(u,'(a)') '#       On ierr /= 0 every T column is written as exactly 0.'
      write(u,'(a)') '# columns:'
      write(u,'(a)') '#  1 lambda[um]  2 a_eff[um]  3 x  4 n  5 k  6 |m-1|  7 ierr  8 t_TM[s]'
      write(u,'(a)') '#  9-11  Qext_R(1:3)   12-14 Qre_R(1:3)   15-17 Qabs_R(1:3)   [Rayleigh]'
      write(u,'(a)') '# 18-20  Qext_T(1:3)   21-23 Qre_T(1:3)   24-26 Qabs_T(1:3)   [T-matrix]'
      write(u,'(a)') '# 27-29  Qext_G(1:3)   30-32 Qre_G(1:3)   33-35 Qabs_G(1:3)   [geom. optics]'
      write(u,'(a)') '# 36 Qext_Mie  37 Qabs_Mie  38 Qsca_R(1)  39 Qsca_T(1)  40 Qsca_G(1)'
   end subroutine write_header


   subroutine forward_amplitude_extinction(wk, a_eff, lambda, m_ref, eps_ba, np, ddelt_in, &
                                           ndgs_in, qext_ori, qre_ori, ierr)
      ! Orientation-resolved extinction and its real-part (birefringence) twin
      ! from the forward-scattering amplitude alone, via the optical theorem
      !     C_ext = (4 pi / k) Im[S_forward],   4 pi / k = 2 lambda,
      ! with the same 4 pi / k prefactor applied to Re[S_forward] to give the
      ! phase-retardation (birefringence) response.  Two AMPL evaluations, no
      ! scattering-sphere integral, so the cost is the T-matrix solve alone.
      ! Same geometry and jori mapping as tmatrix_oriented_cross.
      type(tmatrix_workspace_t), intent(inout) :: wk
      real(wp),    intent(in)  :: a_eff, lambda, eps_ba, ddelt_in
      complex(wp), intent(in)  :: m_ref
      integer,     intent(in)  :: np, ndgs_in
      real(wp),    intent(out) :: qext_ori(3), qre_ori(3)
      integer,     intent(out) :: ierr

      real(wp) :: mrr, mri
      integer  :: ierr_a
      complex(wp) :: vv, vh, hv, hh
      real(wp) :: cext_fac, area
      type(tmatrix_options_t) :: opts
      type(tmatrix_result_t)  :: res
      type(tmatrix_scatmat_t) :: sm

      qext_ori = 0.0_wp;  qre_ori = 0.0_wp
      mrr = real(m_ref, kind=wp)
      mri = abs(aimag(m_ref))
      opts%aspect_ratio = eps_ba
      opts%tolerance    = ddelt_in
      opts%shape        = np
      opts%ndgs         = ndgs_in
      call tmatrix_eval_scattering_matrix(wk, a_eff, lambda, mrr, mri, opts, res, sm)
      if (res%status /= TMATRIX_SUCCESS .and. res%legacy_ierr == 0 .and. &
          res%status /= TMATRIX_ERR_NUMERICAL_NONFINITE) then
         write(error_unit,'(a,a)') ' ERROR: libtmatrix: ', trim(res%message)
         stop 2
      end if
      ierr = res%legacy_ierr
      if (ierr /= 0) return

      cext_fac = 2.0_wp * lambda
      area     = PI * a_eff * a_eff

      call ampl(wk, sm, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, vv, vh, hv, hh, ierr_a)
      if (ierr_a /= 0) then
         write(error_unit,'(a,i0)') ' ERROR: libtmatrix ampl returned ierr = ', ierr_a
         stop 2
      end if
      qext_ori(1) = cext_fac * aimag(vv) / area
      qre_ori(1)  = cext_fac * real(vv, kind=wp) / area
      call ampl(wk, sm, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 90.0_wp, vv, vh, hv, hh, ierr_a)
      if (ierr_a /= 0) then
         write(error_unit,'(a,i0)') ' ERROR: libtmatrix ampl returned ierr = ', ierr_a
         stop 2
      end if
      qext_ori(2) = cext_fac * aimag(vv) / area
      qext_ori(3) = cext_fac * aimag(hh) / area
      qre_ori(2)  = cext_fac * real(vv, kind=wp) / area
      qre_ori(3)  = cext_fac * real(hh, kind=wp) / area
   end subroutine forward_amplitude_extinction


   subroutine read_nodes(fname, lams, as, n)
      character(len=*), intent(in) :: fname
      real(wp), allocatable, intent(out) :: lams(:), as(:)
      integer, intent(out) :: n
      integer :: uu, ios_l, k
      character(len=512) :: buf
      open(newunit=uu, file=fname, status='old', action='read', iostat=ios_l)
      if (ios_l /= 0) then
         write(error_unit,'(a,a)') ' cannot open ', fname
         stop 1
      end if
      n = 0
      do
         read(uu,'(a)',iostat=ios_l) buf
         if (ios_l /= 0) exit
         buf = adjustl(buf)
         if (len_trim(buf) == 0) cycle
         if (buf(1:1) == '#') cycle
         n = n + 1
      end do
      rewind(uu)
      allocate(lams(n), as(n))
      k = 0
      do
         read(uu,'(a)',iostat=ios_l) buf
         if (ios_l /= 0) exit
         buf = adjustl(buf)
         if (len_trim(buf) == 0) cycle
         if (buf(1:1) == '#') cycle
         k = k + 1
         read(buf,*) lams(k), as(k)
      end do
      close(uu)
   end subroutine read_nodes

end program euv_polarized_optics
