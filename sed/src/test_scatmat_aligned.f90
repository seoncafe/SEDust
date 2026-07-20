program test_scatmat_aligned
   !====================================================================
   ! Checks for the aligned-grain polarized-scattering interface
   ! (scatmat_aligned_mod), the reader/query counterpart of calc_polext.
   !
   ! Usage (run from sed/, data paths are relative to ../):
   !   ./test_scatmat_aligned.x [scatmat_table.dat]
   ! With no argument it reads the reduced Stage-2b test table; pass the
   ! production table to validate it.
   !
   ! Anchors:
   !   1. Reader round-trip -- a query at an exact grid node reproduces the
   !      stored file value bit-for-bit (K, F, and Z).
   !   2. Symmetry reconstruction -- the phi mirror and the equatorial mapping
   !      reproduce the sign-flipped stored entry exactly.
   !   3. Eta algebra -- csca_unaligned = Csca_tot - eta Csca_ref exactly, and
   !      the eta=1 unaligned F remainder equals Csca_tot F_tot - Csca_ref F_ref.
   !   4. K consistency -- extinction_matrix_aligned at theta_i = 90, eta = 1
   !      reproduces the f_align-weighted size integrals from the 4-block jori
   !      table (the calc_polext quadrature), and Cpol_ext / Cbir_ext from
   !      dust_extinction agree with K(1,2) / K(3,4).
   !
   ! Each check prints PASS/FAIL with numbers; a FAIL sets a non-zero exit.
   !====================================================================
   use constants,        only: wp, pi, um2cm, deg2rad
   use dust_lib,         only: dust_model_t, build_astrodust, dust_extinction, &
                               dust_set_alignment, dust_set_alignment_profile, &
                               load_scatmat_aligned, free_scatmat_aligned, &
                               scatmat_band, extinction_matrix_aligned, &
                               mueller_matrix_aligned, mueller_matrix_random, &
                               scattering_cross_sections, scm_profile_mismatch, &
                               scm_nband, scm_nti, scm_nts, scm_nphi, scm_ntheta, scm_bytes, &
                               scm_lambda, scm_theta_i, scm_theta_s, scm_phi, scm_theta_ran, &
                               scm_cext_al, scm_cpol_al, scm_cbir_al, scm_csca_al, &
                               scm_csca_tot, scm_csca_ref, &
                               scm_F_tot, scm_F_ref, scm_Z, &
                               scm_profile_name, scm_fmax, scm_a_align, scm_alpha
   use q_table_jori_mod, only: nj_lam, nj_aeff, lam_j, aeff_j, qpol_ext, qbir_ext, &
                               has_bir, falign_hd23
   use size_dist_mod,    only: n_size, a_dist, dn_ad
   implicit none

   character(len=*), parameter :: TEST_TAB = &
      '../tmatrix/output/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400.test.dat'
   character(len=*), parameter :: QTAB = &
      '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: QPOL = &
      '../tmatrix/output/q_astrodust_jori_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: QWAVE = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: QAEFF = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'

   real(wp), parameter :: CM2_TO_UM2 = 1.0e8_wp    ! 1 cm^2 = 1e8 um^2

   character(len=512) :: table
   type(dust_model_t) :: m
   integer :: nfail, st

   nfail = 0
   if (command_argument_count() >= 1) then
      call get_command_argument(1, table)
   else
      table = TEST_TAB
   end if

   write(*,'(a)') '==================================================================='
   write(*,'(a)') ' test_scatmat_aligned'
   write(*,'(a)') '   table = '//trim(table)
   write(*,'(a)') '==================================================================='

   call load_scatmat_aligned(trim(table), st)
   if (st /= 0) then
      write(*,'(a,i0)') ' FATAL: cannot load the scattering table, status = ', st
      stop 2
   end if
   write(*,'(a,i0,a,i0,a,i0,a,i0)') '   bands=', scm_nband, '  nti=', scm_nti, &
        '  nts=', scm_nts, '  nphi=', scm_nphi
   write(*,'(a,f8.2,a)') '   memory measured at load = ', &
        real(scm_bytes,wp)/1.0e6_wp, ' MB'
   write(*,'(a,a,a,f6.3,a,f7.4,a,f5.2)') '   recorded profile: ', &
        trim(scm_profile_name), '  f_max=', scm_fmax, ' a_align=', scm_a_align, &
        ' um alpha=', scm_alpha

   call check_roundtrip(nfail)
   call check_symmetry(nfail)
   call check_eta_algebra(nfail)
   call check_k_consistency(nfail)
   call check_alignment_guard(nfail)

   write(*,'(a)') '-------------------------------------------------------------------'
   if (nfail == 0) then
      write(*,'(a)') ' ALL CHECKS PASSED'
   else
      write(*,'(a,i0,a)') ' ', nfail, ' CHECK(S) FAILED'
   end if
   call free_scatmat_aligned()
   if (nfail /= 0) stop 1

contains

   ! ---- 1. reader round-trip ------------------------------------------
   subroutine check_roundtrip(nf)
      integer, intent(inout) :: nf
      integer  :: ib, ii, is, ip, i, j
      real(wp) :: z(4,4), kmat(4,4), ft(6), fr(6)
      real(wp) :: dz, dk, df, ca, cu
      logical  :: ok

      dz = 0.0_wp;  dk = 0.0_wp;  df = 0.0_wp
      do ib = 1, scm_nband
         ! Z at every node (no folding: all stored nodes are in the fundamental
         ! domain), must reproduce scm_Z bit-for-bit.
         do ii = 1, scm_nti
            do is = 1, scm_nts
               do ip = 1, scm_nphi
                  call mueller_matrix_aligned(ib, scm_theta_i(ii), scm_theta_s(is), &
                                              scm_phi(ip), z)
                  do j = 1, 4
                     do i = 1, 4
                        dz = max(dz, abs(z(i,j) - scm_Z(ii,is,ip,i,j,ib)))
                     end do
                  end do
               end do
            end do
         end do
         ! K at every theta_i node (eta = 1).
         do ii = 1, scm_nti
            call extinction_matrix_aligned(ib, scm_theta_i(ii), 1.0_wp, kmat)
            dk = max(dk, abs(kmat(1,1) - scm_cext_al(ii,ib)))
            dk = max(dk, abs(kmat(1,2) - scm_cpol_al(ii,ib)))
            dk = max(dk, abs(kmat(3,4) - scm_cbir_al(ii,ib)))
            call scattering_cross_sections(ib, scm_theta_i(ii), 1.0_wp, ca, cu)
            dk = max(dk, abs(ca - scm_csca_al(ii,ib)))
         end do
         ! F at every Theta node (the F block's own 1-degree scattering grid).
         do is = 1, scm_ntheta
            call mueller_matrix_random(ib, scm_theta_ran(is), ft, fr)
            do i = 1, 6
               df = max(df, abs(ft(i) - scm_F_tot(is,i,ib)))
               df = max(df, abs(fr(i) - scm_F_ref(is,i,ib)))
            end do
         end do
      end do

      ok = (dz == 0.0_wp .and. dk == 0.0_wp .and. df == 0.0_wp)
      write(*,'(a)')      ' [1] reader round-trip (exact grid nodes reproduce file values)'
      write(*,'(a,es10.2)') '     max |Z - stored|  = ', dz
      write(*,'(a,es10.2)') '     max |K - stored|  = ', dk
      write(*,'(a,es10.2)') '     max |F - stored|  = ', df
      call verdict(ok, nf)
   end subroutine check_roundtrip


   ! ---- 2. symmetry reconstruction ------------------------------------
   subroutine check_symmetry(nf)
      integer, intent(inout) :: nf
      integer  :: ib, ii, is, ip, i, j
      real(wp) :: z_rec(4,4), z_ref(4,4)
      real(wp) :: dphi, dequ
      logical  :: ok

      ib = 1
      dphi = 0.0_wp;  dequ = 0.0_wp
      ! Sample a spread of interior nodes.
      do ii = 1, scm_nti
         do is = 1, scm_nts, max(1, scm_nts/12)
            do ip = 2, scm_nphi-1     ! interior phi so 360-phi is a genuine mirror
               ! (a) phi mirror: query at 360-phi -> stored(phi) with the two
               !     off-diagonal 2x2 blocks flipped.
               call mueller_matrix_aligned(ib, scm_theta_i(ii), scm_theta_s(is), &
                                           360.0_wp - scm_phi(ip), z_rec)
               z_ref = scm_Z(ii,is,ip,:,:,ib)
               call flip_offdiag(z_ref)
               do j = 1, 4
                  do i = 1, 4
                     dphi = max(dphi, abs(z_rec(i,j) - z_ref(i,j)))
                  end do
               end do
            end do
         end do
      end do

      ! (b) equatorial mapping: for theta_i node in (0,90), query
      !     (180-theta_i, 180-theta_s, phi) -> stored(theta_i, theta_s, phi)
      !     flipped. theta_s node reflects onto a theta_s node.
      do ii = 2, scm_nti-1
         do is = 1, scm_nts, max(1, scm_nts/12)
            do ip = 1, scm_nphi
               call mueller_matrix_aligned(ib, 180.0_wp - scm_theta_i(ii), &
                       180.0_wp - scm_theta_s(is), scm_phi(ip), z_rec)
               z_ref = scm_Z(ii,is,ip,:,:,ib)
               call flip_offdiag(z_ref)
               do j = 1, 4
                  do i = 1, 4
                     dequ = max(dequ, abs(z_rec(i,j) - z_ref(i,j)))
                  end do
               end do
            end do
         end do
      end do

      ok = (dphi == 0.0_wp .and. dequ == 0.0_wp)
      write(*,'(a)')      ' [2] symmetry reconstruction (reconstructed == sign-flipped stored)'
      write(*,'(a,es10.2)') '     phi mirror   max |diff| = ', dphi
      write(*,'(a,es10.2)') '     equatorial   max |diff| = ', dequ
      call verdict(ok, nf)
   end subroutine check_symmetry


   ! ---- 3. eta algebra ------------------------------------------------
   subroutine check_eta_algebra(nf)
      integer, intent(inout) :: nf
      integer  :: ib, is, i, k
      real(wp) :: eta(2), e, ca, cu, resid, maxresid
      real(wp) :: fu_routine, fu_direct, dfu, ft(6), fr(6)
      real(wp) :: closure, half_int
      logical  :: ok

      eta = [1.0_wp, 0.37_wp]
      ib = 1
      maxresid = 0.0_wp;  dfu = 0.0_wp
      do k = 1, 2
         e = eta(k)
         ! Scalar: csca_unaligned = Csca_tot - eta Csca_ref, so eta Csca_ref +
         ! csca_unaligned must return Csca_tot to rounding. (csca_aligned is the
         ! direction-resolved eta Csca_al(theta_i); it recovers Csca_tot only
         ! after averaging over incidence, so the EXACT identity is on the band
         ! scalar Csca_ref.)
         call scattering_cross_sections(ib, 45.0_wp, e, ca, cu)
         resid = abs((e*scm_csca_ref(ib) + cu) - scm_csca_tot(ib)) / scm_csca_tot(ib)
         maxresid = max(maxresid, resid)
         if (abs(cu - (scm_csca_tot(ib) - e*scm_csca_ref(ib))) /= 0.0_wp) &
            maxresid = max(maxresid, 1.0_wp)   ! subtraction must be bit-exact

         ! F remainder: F_unal(eta) = Csca_tot F_tot - eta Csca_ref F_ref must
         ! match the direct combination of the stored matrices at every node.
         do is = 1, scm_ntheta
            call mueller_matrix_random(ib, scm_theta_ran(is), ft, fr)
            do i = 1, 6
               fu_routine = scm_csca_tot(ib)*ft(i) - e*scm_csca_ref(ib)*fr(i)
               fu_direct  = scm_csca_tot(ib)*scm_F_tot(is,i,ib) &
                          - e*scm_csca_ref(ib)*scm_F_ref(is,i,ib)
               dfu = max(dfu, abs(fu_routine - fu_direct))
            end do
         end do
      end do

      ! Physical closure of the eta = 1 unaligned remainder: (1/2) INT F_unal_11
      ! dcos over the F Theta grid vs Csca_tot - Csca_ref.
      half_int = 0.0_wp
      do is = 2, scm_ntheta
         half_int = half_int + trap_cos(ib, is)
      end do
      half_int = 0.5_wp * half_int
      closure = half_int / (scm_csca_tot(ib) - scm_csca_ref(ib))

      ok = (maxresid <= 1.0e-12_wp .and. dfu == 0.0_wp)
      write(*,'(a)')      ' [3] eta algebra (eta = 1, 0.37)'
      write(*,'(a,es10.2)') '     max rel |eta Csca_ref + csca_unal - Csca_tot| = ', maxresid
      write(*,'(a,es10.2)') '     max |F_unal(routine) - F_unal(direct)|        = ', dfu
      write(*,'(a,f9.5)')   '     eta=1 remainder closure (1/2 INT F_unal_11 dcos)/(Csca_tot-Csca_ref) = ', closure
      call verdict(ok, nf)
   end subroutine check_eta_algebra


   real(wp) function trap_cos(ib, is) result(seg)
      ! One trapezoid segment of INT F_unal_11 dcos(Theta) between Theta nodes
      ! is-1 and is, at eta = 1 (F_unal_11 = Csca_tot F_tot_11 - Csca_ref F_ref_11).
      integer, intent(in) :: ib, is
      real(wp) :: c0, c1, f0, f1
      c0 = cos(scm_theta_ran(is-1)*deg2rad)
      c1 = cos(scm_theta_ran(is)  *deg2rad)
      f0 = scm_csca_tot(ib)*scm_F_tot(is-1,1,ib) - scm_csca_ref(ib)*scm_F_ref(is-1,1,ib)
      f1 = scm_csca_tot(ib)*scm_F_tot(is,  1,ib) - scm_csca_ref(ib)*scm_F_ref(is,  1,ib)
      seg = 0.5_wp * (f0 + f1) * (c0 - c1)      ! dcos > 0 as Theta increases
   end function trap_cos


   ! ---- 4. K consistency ----------------------------------------------
   subroutine check_k_consistency(nf)
      integer, intent(inout) :: nf
      integer  :: ib, nlam_m
      logical  :: exact
      real(wp) :: lam_band, kmat(4,4)
      real(wp) :: kpol_um2, kbir_um2               ! K(1,2), K(3,4) at 90 deg, um^2/H
      real(wp) :: cpol_int, cbir_int               ! calc_polext-style integrals, cm^2/H
      real(wp) :: rel_pol_jori, rel_bir_jori
      real(wp) :: rel_pol_de, rel_bir_de
      real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), Cpol_ext(:), Cbir_ext(:), lam_m(:)
      logical  :: ok

      ! Build the model with the 4-block jori table so dust_extinction returns a
      ! non-zero Cbir_ext; this also leaves qpol_ext / qbir_ext / lam_j / aeff_j
      ! and the size distribution loaded for the independent quadrature.
      call build_astrodust(m, QTAB, SIZED, 100, 2.7_wp, 5.0e3_wp, status=st, &
                           qpol_path=QPOL, qpol_wave_path=QWAVE, qpol_aeff_path=QAEFF)
      if (st /= 0) then
         write(*,'(a,i0)') ' [4] K consistency: build_astrodust failed, status = ', st
         nf = nf + 1;  return
      end if
      if (.not. has_bir) then
         write(*,'(a)') ' [4] K consistency: jori table has no birefringence block (unexpected)'
         nf = nf + 1;  return
      end if

      nlam_m = size(m%lam)
      allocate(lam_m(nlam_m), Cext(nlam_m), Cabs(nlam_m), Csca(nlam_m), &
               Cpol_ext(nlam_m), Cbir_ext(nlam_m))
      lam_m = m%lam
      call dust_extinction(m, Cext, Cabs, Csca, Cpol_ext=Cpol_ext, Cbir_ext=Cbir_ext)

      ! Pick a band; use 0.55 um, which is one of the shipped bands.
      call scatmat_band(0.55_wp, ib, exact)
      lam_band = scm_lambda(ib)

      ! K(1,2), K(3,4) at theta_i = 90, eta = 1, converted to cm^2/H below.
      call extinction_matrix_aligned(ib, 90.0_wp, 1.0_wp, kmat)
      kpol_um2 = kmat(1,2)
      kbir_um2 = kmat(3,4)

      ! Independent f_align-weighted size integral from the jori Q table,
      ! interpolated in log-wavelength to exactly lam_band -- the band is a UBVRI
      ! effective wavelength, not a DH21 grid node, and the generator computed K
      ! at lam_band, so a nearest-node integral would carry a spurious
      ! wavelength offset (the calc_polext quadrature; cm^2/H).
      cpol_int = polext_at_lambda(lam_band, qpol_ext)
      cbir_int = polext_at_lambda(lam_band, qbir_ext)

      ! dust_extinction interpolated in log-wavelength to lam_band (cm^2/H).
      rel_pol_jori = reldiff(cpol_int*CM2_TO_UM2, kpol_um2)
      rel_bir_jori = reldiff(cbir_int*CM2_TO_UM2, kbir_um2)
      rel_pol_de   = reldiff(interp_lam(lam_m, nlam_m, Cpol_ext, lam_band)*CM2_TO_UM2, kpol_um2)
      rel_bir_de   = reldiff(interp_lam(lam_m, nlam_m, Cbir_ext, lam_band)*CM2_TO_UM2, kbir_um2)

      write(*,'(a)')        ' [4] K consistency at theta_i = 90, eta = 1'
      write(*,'(a,f6.3,a,l1,a)') '     band lambda = ', lam_band, ' um (exact match = ', exact, ')'
      write(*,'(a,es12.4,a,es12.4)') '     K(1,2)=Cpol_al [um^2/H] = ', kpol_um2, &
           '   K(3,4)=Cbir_al [um^2/H] = ', kbir_um2
      write(*,'(a,es12.4,a,es12.4)') '     jori integral  Cpol     = ', cpol_int*CM2_TO_UM2, &
           '   Cbir     = ', cbir_int*CM2_TO_UM2
      write(*,'(a,es12.4,a,es12.4)') '     dust_extinction Cpol_ext = ', &
           interp_lam(lam_m, nlam_m, Cpol_ext, lam_band)*CM2_TO_UM2, &
           '   Cbir_ext = ', interp_lam(lam_m, nlam_m, Cbir_ext, lam_band)*CM2_TO_UM2
      write(*,'(a,es10.2,a,es10.2)') '     rel diff jori vs K:  Cpol = ', rel_pol_jori, &
           '   Cbir = ', rel_bir_jori
      write(*,'(a,es10.2,a,es10.2)') '     rel diff d_ext vs K: Cpol = ', rel_pol_de, &
           '   Cbir = ', rel_bir_de

      ! dust_extinction and the jori quadrature share the size grid and f_align,
      ! so they must agree with K(1,2)/K(3,4) tightly (the ~1e-3 target for the
      ! dust_extinction comparison); the jori quadrature is an independent
      ! integral and is allowed slightly more slack.
      ok = (rel_pol_de <= 1.0e-3_wp .and. rel_bir_de <= 1.0e-3_wp .and. &
            rel_pol_jori <= 1.0e-2_wp .and. rel_bir_jori <= 1.0e-2_wp)
      call verdict(ok, nf)

      deallocate(lam_m, Cext, Cabs, Csca, Cpol_ext, Cbir_ext)
   end subroutine check_k_consistency


   ! ---- 5. alignment-consistency guard --------------------------------
   subroutine check_alignment_guard(nf)
      ! With the aligned scattering table loaded, dust_set_alignment must accept
      ! the recorded profile silently (status 0) and flag any departure (status 4
      ! + scm_profile_mismatch), and dust_set_alignment_profile must flag a
      ! tabulated profile. Uses the model built in check_k_consistency.
      integer, intent(inout) :: nf
      integer  :: s_same, s_diff, s_tab
      real(wp) :: aeff_in(3), fal_in(3)
      logical  :: ok

      ! (a) same profile as the table -> accepted.
      call dust_set_alignment(m, scm_fmax, scm_a_align, scm_alpha, status=s_same)
      ok = (s_same == 0 .and. .not. scm_profile_mismatch)

      ! (b) a different f_max -> flagged, non-fatal.
      call dust_set_alignment(m, 0.5_wp*scm_fmax, scm_a_align, scm_alpha, status=s_diff)
      ok = ok .and. (s_diff == 4 .and. scm_profile_mismatch)

      ! (c) a tabulated profile -> flagged.
      aeff_in = [0.01_wp, 0.1_wp, 1.0_wp]
      fal_in  = [0.0_wp, 0.5_wp, 1.0_wp]
      call dust_set_alignment_profile(m, aeff_in, fal_in, status=s_tab)
      ok = ok .and. (s_tab == 4 .and. scm_profile_mismatch)

      ! restore the recorded profile for hygiene
      call dust_set_alignment(m, scm_fmax, scm_a_align, scm_alpha)

      write(*,'(a)')       ' [5] alignment-consistency guard'
      write(*,'(a,i0,a,i0,a,i0)') '     status: same profile = ', s_same, &
           '   changed f_max = ', s_diff, '   tabulated = ', s_tab
      call verdict(ok, nf)
   end subroutine check_alignment_guard


   real(wp) function polext_at_lambda(lam0, qarr) result(cint)
      ! calc_polext-style size integral at exactly lam0: the integral is formed
      ! at the two bracketing jori wavelength nodes and interpolated log-linearly
      ! in wavelength, matching the way the generator sampled the band [cm^2/H].
      real(wp), intent(in) :: lam0
      real(wp), intent(in) :: qarr(:,:)      ! (nj_lam, nj_aeff)
      integer  :: jlo, jhi
      real(wp) :: t, clo, chi
      call ibracket_log(lam_j, nj_lam, lam0, jlo, t)
      jhi = jlo + 1
      clo = polext_integral(jlo, qarr)
      chi = polext_integral(jhi, qarr)
      cint = (1.0_wp - t)*clo + t*chi
   end function polext_at_lambda


   real(wp) function polext_integral(jw, qarr) result(cint)
      ! sum_a dn(a) f_align(a) Q(a) pi a^2, with Q interpolated log-linearly in
      ! a onto the size-distribution grid at fixed wavelength index jw [cm^2/H].
      integer,  intent(in) :: jw
      real(wp), intent(in) :: qarr(:,:)      ! (nj_lam, nj_aeff)
      integer  :: ia
      real(wp) :: qi
      cint = 0.0_wp
      do ia = 1, n_size
         qi = interp_a(jw, a_dist(ia), qarr)
         cint = cint + dn_ad(ia) * falign_hd23(a_dist(ia)) * qi &
                     * pi * (a_dist(ia)*um2cm)**2
      end do
   end function polext_integral


   real(wp) function interp_lam(g, n, y, x) result(v)
      ! Linear-in-log-wavelength interpolation of y(g) at x, clamped at the ends.
      real(wp), intent(in) :: g(:), y(:)
      integer,  intent(in) :: n
      real(wp), intent(in) :: x
      integer  :: jlo
      real(wp) :: t
      call ibracket_log(g, n, x, jlo, t)
      v = (1.0_wp - t)*y(jlo) + t*y(jlo+1)
   end function interp_lam


   subroutine ibracket_log(g, n, x, jlo, t)
      ! Bracket x in the strictly increasing grid g(1:n); linear weight in log g.
      real(wp), intent(in)  :: g(:)
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: x
      integer,  intent(out) :: jlo
      real(wp), intent(out) :: t
      integer :: lo, hi, mid
      if (x <= g(1)) then
         jlo = 1;  t = 0.0_wp;  return
      end if
      if (x >= g(n)) then
         jlo = n-1;  t = 1.0_wp;  return
      end if
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (g(mid) <= x) then;  lo = mid;  else;  hi = mid;  end if
      end do
      jlo = lo
      t = (log(x) - log(g(jlo))) / (log(g(jlo+1)) - log(g(jlo)))
   end subroutine ibracket_log


   real(wp) function interp_a(jw, a_target, qarr) result(q)
      integer,  intent(in) :: jw
      real(wp), intent(in) :: a_target
      real(wp), intent(in) :: qarr(:,:)
      integer  :: l, h, mm
      real(wp) :: xa, ta
      if (a_target <= aeff_j(1)) then
         q = qarr(jw, 1);        return
      end if
      if (a_target >= aeff_j(nj_aeff)) then
         q = qarr(jw, nj_aeff);  return
      end if
      xa = log(a_target)
      l = 1;  h = nj_aeff
      do while (h - l > 1)
         mm = (l + h) / 2
         if (log(aeff_j(mm)) <= xa) then
            l = mm
         else
            h = mm
         end if
      end do
      ta = (xa - log(aeff_j(l))) / (log(aeff_j(h)) - log(aeff_j(l)))
      q  = (1.0_wp - ta) * qarr(jw, l) + ta * qarr(jw, h)
   end function interp_a


   ! ---- small utilities -----------------------------------------------
   real(wp) function reldiff(a, b) result(r)
      real(wp), intent(in) :: a, b
      real(wp) :: s
      s = max(abs(a), abs(b), tiny(1.0_wp))
      r = abs(a - b) / s
   end function reldiff

   subroutine flip_offdiag(z)
      real(wp), intent(inout) :: z(4,4)
      z(1,3) = -z(1,3);  z(1,4) = -z(1,4)
      z(2,3) = -z(2,3);  z(2,4) = -z(2,4)
      z(3,1) = -z(3,1);  z(3,2) = -z(3,2)
      z(4,1) = -z(4,1);  z(4,2) = -z(4,2)
   end subroutine flip_offdiag

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

end program test_scatmat_aligned
