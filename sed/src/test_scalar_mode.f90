program test_scalar_mode
   !====================================================================
   ! Checks for the explicit scalar-only astrodust build mode
   ! (load_polarized_optics = .false.) and the polarized-optics query
   ! dust_has_polarized_optics.
   !
   ! Usage (run from sed/, data paths are relative to ../):
   !   ./test_scalar_mode.x
   !
   ! Anchors:
   !   1. Scalar vs polarized extinction -- Cext, Cabs, Csca, gbar over the
   !      full wavelength grid are identical (scattering optics are attached in
   !      both builds; only the polarized channel differs).
   !   2. Scalar vs polarized SED -- the total lambda*I_lambda for one Mathis
   !      field is identical (the total emission never touches the polarized
   !      optics, and the default heuristic solver is serial/deterministic).
   !   3. Scalar polarized outputs -- Cpol_ext and Cbir_ext are exactly zero,
   !      and dust_has_polarized_optics is .false. for the scalar model and
   !      .true. for the polarized one.
   !   4. No decompression scratch file -- a scalar build opens no polarized Q
   !      table, so no q_jori_<pid>.dat is created in TMPDIR/CWD; the polarized
   !      build creates and then removes one, so none lingers either.
   !   5. Contradiction -- load_polarized_optics = .false. combined with an
   !      explicit qpol_path (or scatmat_path) returns the documented status 5.
   !
   ! Each check prints PASS/FAIL with numbers; a FAIL sets a non-zero exit.
   !====================================================================
   use constants, only: wp
   use radfield,  only: J_Mathis
   use dust_lib,  only: dust_model_t, build_astrodust, dust_extinction, &
                        dust_emission, dust_has_polarized_optics, &
                        dust_set_alignment, dust_nlam
   implicit none

   character(len=*), parameter :: QTAB  = &
      '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'
   ! The default polarized table, passed only to exercise the contradiction
   ! path (its contents are never read there -- the clash is caught first).
   character(len=*), parameter :: QPOL  = &
      '../data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz'

   integer,  parameter :: NT_IN = 100
   real(wp), parameter :: T_LO = 2.7_wp, T_HI = 5.0e3_wp
   ! HD23 Table 1 alignment parameters, used to force the polarized build's
   ! alignment to zero (f_max = 0) so the two builds are maximally comparable.
   real(wp), parameter :: A_ALIGN = 0.0749_wp, ALPHA_ALIGN = 1.80_wp

   ! The compared quantities are algebraically identical between the two builds,
   ! so the measured differences are 0; the tolerances are a regression guard.
   real(wp), parameter :: TOL_EXT = 1.0e-14_wp
   real(wp), parameter :: TOL_SED = 1.0e-14_wp

   type(dust_model_t) :: m_pol, m_scalar, m_tmp, m_bad
   integer :: nfail, st, nlam
   real(wp), allocatable :: Cext_p(:), Cabs_p(:), Csca_p(:), gbar_p(:)
   real(wp), allocatable :: Cext_s(:), Cabs_s(:), Csca_s(:), gbar_s(:)
   real(wp), allocatable :: Cpol_ext_s(:), Cbir_ext_s(:)
   real(wp), allocatable :: J(:), total_p(:), total_s(:)
   logical :: has_pol_p, has_pol_s
   character(len=512) :: tmpdir
   integer :: envstat

   nfail = 0
   call get_environment_variable('TMPDIR', tmpdir, status=envstat)
   if (envstat /= 0 .or. len_trim(tmpdir) == 0) tmpdir = '/tmp'

   write(*,'(a)') '==================================================================='
   write(*,'(a)') ' test_scalar_mode'
   write(*,'(a)') '   qtable = '//QTAB
   write(*,'(a)') '   TMPDIR = '//trim(tmpdir)
   write(*,'(a)') '==================================================================='

   ! ---- polarized build (default: loads the polarized optics) ----------
   call build_astrodust(m_pol, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: polarized build_astrodust failed, status = ', st
      stop 2
   end if
   ! Force alignment to zero so nothing but the presence of polarized optics
   ! distinguishes the two builds (alignment enters neither extinction scalars
   ! nor the total SED; this just removes any confounder).
   call dust_set_alignment(m_pol, 0.0_wp, A_ALIGN, ALPHA_ALIGN, status=st)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: dust_set_alignment failed, status = ', st
      stop 2
   end if
   has_pol_p = dust_has_polarized_optics(m_pol)

   nlam = dust_nlam(m_pol)
   allocate(Cext_p(nlam), Cabs_p(nlam), Csca_p(nlam), gbar_p(nlam))
   allocate(Cext_s(nlam), Cabs_s(nlam), Csca_s(nlam), gbar_s(nlam))
   allocate(Cpol_ext_s(nlam), Cbir_ext_s(nlam))
   allocate(J(nlam), total_p(nlam), total_s(nlam))

   call dust_extinction(m_pol, Cext_p, Cabs_p, Csca_p, gbar=gbar_p)
   ! One cell's field: the Mathis ISRF at U = 1.
   call J_Mathis(1.0_wp, m_pol%lam, J)
   call dust_emission(m_pol, J, total_p)

   ! ---- scalar-only build (never opens the polarized table) ------------
   call build_astrodust(m_scalar, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st, &
                        load_polarized_optics=.false.)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: scalar build_astrodust failed, status = ', st
      stop 2
   end if
   has_pol_s = dust_has_polarized_optics(m_scalar)

   call dust_extinction(m_scalar, Cext_s, Cabs_s, Csca_s, gbar=gbar_s, &
                        Cpol_ext=Cpol_ext_s, Cbir_ext=Cbir_ext_s)
   call dust_emission(m_scalar, J, total_s)

   call check_extinction(nfail)
   call check_sed(nfail)
   call check_polarized_outputs(nfail)
   call check_no_scratch(nfail)
   call check_contradiction(nfail)

   write(*,'(a)') '-------------------------------------------------------------------'
   if (nfail == 0) then
      write(*,'(a)') ' ALL CHECKS PASSED'
   else
      write(*,'(a,i0,a)') ' ', nfail, ' CHECK(S) FAILED'
   end if
   if (nfail /= 0) stop 1

contains

   ! ---- 1. extinction scalars identical ------------------------------
   subroutine check_extinction(nf)
      integer, intent(inout) :: nf
      real(wp) :: dext, dabs, dsca, dg
      logical  :: ok
      dext = maxreldiff(Cext_p, Cext_s)
      dabs = maxreldiff(Cabs_p, Cabs_s)
      dsca = maxreldiff(Csca_p, Csca_s)
      dg   = maxreldiff(gbar_p, gbar_s)
      ok = (dext <= TOL_EXT .and. dabs <= TOL_EXT .and. &
            dsca <= TOL_EXT .and. dg <= TOL_EXT)
      write(*,'(a)')        ' [1] extinction scalars identical (scalar vs polarized build)'
      write(*,'(a,es10.2,a,es10.2)') '     max rel |dCext| = ', dext, '   |dCabs| = ', dabs
      write(*,'(a,es10.2,a,es10.2)') '     max rel |dCsca| = ', dsca, '   |dgbar| = ', dg
      write(*,'(a,es10.2)') '     tol = ', TOL_EXT
      call verdict(ok, nf)
   end subroutine check_extinction

   ! ---- 2. total SED identical ---------------------------------------
   subroutine check_sed(nf)
      integer, intent(inout) :: nf
      real(wp) :: dsed
      logical  :: ok
      dsed = maxreldiff(total_p, total_s)
      ok = (dsed <= TOL_SED)
      write(*,'(a)')        ' [2] total SED identical (Mathis U=1; heuristic solver, serial)'
      write(*,'(a,es10.2,a,es10.2)') '     max rel |d(lamI)| = ', dsed, '   tol = ', TOL_SED
      call verdict(ok, nf)
   end subroutine check_sed

   ! ---- 3. scalar polarized outputs zero + query -----------------------
   subroutine check_polarized_outputs(nf)
      integer, intent(inout) :: nf
      real(wp) :: mpol, mbir
      logical  :: ok
      mpol = maxval(abs(Cpol_ext_s))
      mbir = maxval(abs(Cbir_ext_s))
      ok = (mpol == 0.0_wp .and. mbir == 0.0_wp .and. &
            has_pol_p .and. (.not. has_pol_s))
      write(*,'(a)')        ' [3] scalar polarized outputs and dust_has_polarized_optics'
      write(*,'(a,es10.2,a,es10.2)') '     scalar max|Cpol_ext| = ', mpol, &
           '   max|Cbir_ext| = ', mbir
      write(*,'(a,l1,a,l1)') '     has_polarized: polarized build = ', has_pol_p, &
           '   scalar build = ', has_pol_s
      call verdict(ok, nf)
   end subroutine check_polarized_outputs

   ! ---- 4. no decompression scratch file for a scalar build -----------
   subroutine check_no_scratch(nf)
      integer, intent(inout) :: nf
      logical :: before, after
      logical :: ok
      integer :: st4
      ! The polarized build above created and then removed a q_jori_<pid>.dat,
      ! so nothing should linger now.
      before = scratch_present(tmpdir)
      ! A fresh scalar build must open no polarized table, hence create no
      ! scratch file at all.
      call build_astrodust(m_tmp, QTAB, SIZED, NT_IN, T_LO, T_HI, status=st4, &
                           load_polarized_optics=.false.)
      after = scratch_present(tmpdir)
      ok = (st4 == 0 .and. (.not. before) .and. (.not. after))
      write(*,'(a)')        ' [4] no scratch file q_jori_*.dat for a scalar build'
      write(*,'(a,l1,a,l1,a,i0)') '     lingering before scalar build = ', before, &
           '   after = ', after, '   build status = ', st4
      call verdict(ok, nf)
   end subroutine check_no_scratch

   ! ---- 5. contradiction: scalar mode + explicit polarized path -------
   subroutine check_contradiction(nf)
      integer, intent(inout) :: nf
      integer :: s_qpol, s_scat
      logical :: ok
      ! (a) explicit qpol_path with load_polarized_optics = .false.
      call build_astrodust(m_bad, QTAB, SIZED, NT_IN, T_LO, T_HI, status=s_qpol, &
                           qpol_path=QPOL, load_polarized_optics=.false.)
      ! (b) explicit scatmat_path with load_polarized_optics = .false.
      call build_astrodust(m_bad, QTAB, SIZED, NT_IN, T_LO, T_HI, status=s_scat, &
                           scatmat_path=QPOL, load_polarized_optics=.false.)
      ok = (s_qpol == 5 .and. s_scat == 5)
      write(*,'(a)')        ' [5] contradiction (scalar mode + explicit polarized path)'
      write(*,'(a,i0,a,i0)') '     status: qpol_path = ', s_qpol, &
           '   scatmat_path = ', s_scat
      call verdict(ok, nf)
   end subroutine check_contradiction

   ! ---- utilities -----------------------------------------------------
   real(wp) function maxreldiff(a, b) result(d)
      ! Maximum over the grid of |a-b| / max(|a|,|b|,tiny). Exactly 0 when the
      ! two arrays are bit-identical (the expected result here).
      real(wp), intent(in) :: a(:), b(:)
      integer  :: k
      real(wp) :: s
      d = 0.0_wp
      do k = 1, size(a)
         s = max(abs(a(k)), abs(b(k)), tiny(1.0_wp))
         d = max(d, abs(a(k) - b(k)) / s)
      end do
   end function maxreldiff

   logical function scratch_present(dir) result(present_any)
      ! .true. iff any q_jori_*.dat exists in dir or the current directory.
      character(len=*), intent(in) :: dir
      integer :: es, cs
      call execute_command_line('ls '//trim(dir)//'/q_jori_*.dat ./q_jori_*.dat '// &
                                '>/dev/null 2>&1', exitstat=es, cmdstat=cs, wait=.true.)
      present_any = (cs == 0 .and. es == 0)
   end function scratch_present

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

end program test_scalar_mode
