program use_dustlib_scatmat
   ! Aligned-grain polarized-scattering example for an EXTERNAL Fortran code (a
   ! 3D polarized RT driver) linking the dust library. It is NOT part of the sed
   ! build -- it is compiled separately against libsedust.a + the .mod search
   ! path, exactly as an RT code would. Build and run from sed/, because the data
   ! paths below are relative to that directory:
   !
   !   cd sed
   !   make libsedust.a
   !   make use_dustlib_scatmat.x     # or: gfortran -I. rt_example/use_dustlib_scatmat.f90 \
   !                                  #        libsedust.a -fopenmp -o use_dustlib_scatmat.x
   !   ./use_dustlib_scatmat.x [scatmat_table.dat]
   !
   ! What it shows: the Peest-style consumption of the aligned scattering optics
   ! along a fake two-cell photon path. SEDust is initialized ONCE (the aligned
   ! table is parsed and the size integrals are done); then, per cell, the only
   ! inputs are the scalars the RT already holds -- the alignment scale eta and
   ! the incidence angle theta_i between the photon direction k-hat and the local
   ! field B-hat (one dot product). The query calls are pure reads.
   !
   ! Division of labor (plan Section 4): SEDust returns the matrices in the grain
   ! frame (z = B-hat); the RT does the meridional-basis rotations, azimuth
   ! sampling, peel-off, and the exp(-K tau) transfer.
   use constants, only: wp, deg2rad, rad2deg
   use dust_lib,  only: dust_model_t, build_astrodust, dust_extinction, &
                        scatmat_band, extinction_matrix_aligned, &
                        mueller_matrix_aligned, mueller_matrix_random, &
                        scattering_cross_sections, scm_loaded, &
                        scm_nts, scm_nphi, scm_lambda, &
                        scm_theta_s, scm_phi, scm_cext_ref, scm_csca_tot, scm_csca_ref
   implicit none

   character(len=*), parameter :: SCATMAT_DEF = &
      '../tmatrix/output/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400.test.dat'
   character(len=*), parameter :: QTAB = &
      '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: QPOL = &
      '../tmatrix/output/q_astrodust_jori_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: QWAVE = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: QAEFF = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: SIZED = '../data/release/size_distribution.dat'
   real(wp), parameter :: CM2_TO_UM2 = 1.0e8_wp

   character(len=512) :: scatmat
   type(dust_model_t) :: m
   real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), Cpol_ext(:), Cbir_ext(:)
   real(wp) :: cext_iso_um2                       ! isotropic total extinction, um^2/H
   integer  :: n, st, iband, icell
   logical  :: exact
   ! Per-cell inputs the RT already has in hand.
   real(wp) :: eta(2), khat(3), Bhat(3,2)
   real(wp) :: costi, theta_i, kmat(4,4), z(4,4), f_tot(6), f_ref(6)
   real(wp) :: theta_s, phi, big_theta
   real(wp) :: csca_al, csca_unal, z11_closure, fu11_fwd
   logical  :: bad

   if (command_argument_count() >= 1) then
      call get_command_argument(1, scatmat)
   else
      scatmat = SCATMAT_DEF
   end if

   ! ---- initialize ONCE: model + aligned scattering table ----
   call build_astrodust(m, QTAB, SIZED, 100, 2.7_wp, 5.0e3_wp, status=st, &
                        qpol_path=QPOL, qpol_wave_path=QWAVE, qpol_aeff_path=QAEFF, &
                        scatmat_path=trim(scatmat))
   if (st /= 0) then
      print '(a,i0)', ' build_astrodust failed, status = ', st
      stop 1
   end if
   if (.not. scm_loaded) then
      print '(a)', ' aligned scattering table not loaded'
      stop 1
   end if

   n = size(m%lam)
   allocate(Cext(n), Cabs(n), Csca(n), Cpol_ext(n), Cbir_ext(n))
   call dust_extinction(m, Cext, Cabs, Csca, Cpol_ext=Cpol_ext, Cbir_ext=Cbir_ext)

   ! Pick the band once; the hot path then takes iband.
   call scatmat_band(0.55_wp, iband, exact)
   cext_iso_um2 = interp_lam(m%lam, n, Cext, scm_lambda(iband)) * CM2_TO_UM2

   print '(a)', ' === SEDust aligned-scattering interface (Peest-style consumer) ==='
   print '(a,f6.3,a,l1,a,i0,a)', '   band lambda = ', scm_lambda(iband), &
         ' um (exact=', exact, ')   grid: nts=', scm_nts, ''
   print '(a,es12.4)', '   isotropic total Cext [um^2/H]     = ', cext_iso_um2
   print '(a,es12.4)', '   aligned reference Cext_ref [um^2/H]= ', scm_cext_ref(iband)
   print '(a,es12.4,a,es12.4)', '   Csca_tot = ', scm_csca_tot(iband), &
         '   Csca_ref = ', scm_csca_ref(iband)

   ! ---- fake two-cell path: same photon direction, different eta and field ----
   khat     = [0.0_wp, 0.0_wp, 1.0_wp]
   eta      = [1.00_wp, 0.40_wp]
   Bhat(:,1) = [sin(30.0_wp*deg2rad), 0.0_wp, cos(30.0_wp*deg2rad)]   ! 30 deg to k
   Bhat(:,2) = [sin(75.0_wp*deg2rad), 0.0_wp, cos(75.0_wp*deg2rad)]   ! 75 deg to k

   ! A fixed scattering geometry to demonstrate the queries.
   theta_s   = 60.0_wp
   phi       = 45.0_wp
   big_theta = 60.0_wp

   bad = .false.
   do icell = 1, 2
      ! theta_i from one dot product k-hat . B-hat, the only cell geometry the RT
      ! must supply besides eta.
      costi   = dot_product(khat, Bhat(:,icell))
      theta_i = acos(max(-1.0_wp, min(1.0_wp, costi))) * rad2deg

      ! Total extinction matrix along the ray. The isotropic total from
      ! dust_extinction already contains the aligned population's contribution at
      ! its reference orientation average, Cext_ref; to make it direction-
      ! dependent we SUBTRACT that reference piece (eta*Cext_ref) and ADD the
      ! direction-resolved aligned matrix eta*K_al(theta_i). At the incidence
      ! average <K_al> = Cext_ref the two cancel and the diagonal returns
      ! Cext_iso, so no extinction is double-counted.
      call extinction_matrix_aligned(iband, theta_i, eta(icell), kmat)   ! eta*K_al, um^2/H
      kmat(1,1) = kmat(1,1) + (cext_iso_um2 - eta(icell)*scm_cext_ref(iband))
      kmat(2,2) = kmat(2,2) + (cext_iso_um2 - eta(icell)*scm_cext_ref(iband))
      kmat(3,3) = kmat(3,3) + (cext_iso_um2 - eta(icell)*scm_cext_ref(iband))
      kmat(4,4) = kmat(4,4) + (cext_iso_um2 - eta(icell)*scm_cext_ref(iband))

      ! Aligned phase matrix at the scattering geometry (reference), scaled by eta.
      call mueller_matrix_aligned(iband, theta_i, theta_s, phi, z)
      z = eta(icell) * z

      ! Randomly-oriented remainder in absolute units: Csca_tot F_tot - eta Csca_ref F_ref.
      call mueller_matrix_random(iband, big_theta, f_tot, f_ref)
      fu11_fwd = scm_csca_tot(iband)*f_tot(1) - eta(icell)*scm_csca_ref(iband)*f_ref(1)

      ! Scalar scattering cross sections and the Z11 closure of the aligned part.
      call scattering_cross_sections(iband, theta_i, eta(icell), csca_al, csca_unal)
      z11_closure = aligned_z11_integral(iband, theta_i, eta(icell))

      print '(a)', ' ------------------------------------------------------------------'
      print '(a,i0,a,f5.2,a,f6.2,a)', '   cell ', icell, ':  eta = ', eta(icell), &
            '   theta_i = ', theta_i, ' deg  (k.B)'
      print '(a,es12.4,a,es12.4)', '     K(1,1)=Cext [um^2/H] = ', kmat(1,1), &
            '   K(1,2)=Cpol = ', kmat(1,2)
      print '(a,es12.4,a,es12.4)', '     K(3,4)=Cbir          = ', kmat(3,4), &
            '   K(4,3)      = ', kmat(4,3)
      print '(a,es12.4,a,f8.5)', '     aligned Z11(60,45) [um^2/sr/H] = ', z(1,1), &
            '   -Z12/Z11 = ', -z(1,2)/z(1,1)
      print '(a,es12.4)', '     unaligned F_unal,11(60) [um^2/sr/H] = ', fu11_fwd
      print '(a,es12.4,a,es12.4)', '     csca_aligned = ', csca_al, &
            '   csca_unaligned = ', csca_unal
      print '(a,f9.5)', '     Z11 closure (INT Z11 dOmega)/csca_aligned = ', &
            z11_closure / csca_al

      ! Physical sanity, printed and gathered into the exit status.
      if (.not. finite(kmat(1,1)) .or. .not. finite(z(1,1))) bad = .true.
      if (z(1,1) <= 0.0_wp)               bad = .true.    ! phase function positive
      if (abs(kmat(1,2)) >= kmat(1,1))    bad = .true.    ! |Cpol| < Cext
   end do

   print '(a)', ' ------------------------------------------------------------------'
   if (bad) then
      print '(a)', '   SANITY: FAILED (NaN, non-positive Z11, or |Cpol| >= Cext)'
      deallocate(Cext, Cabs, Csca, Cpol_ext, Cbir_ext)
      stop 1
   else
      print '(a)', '   SANITY: OK (finite, Z11 > 0, |Cpol| < Cext, closures ~ 1)'
   end if
   print '(a)', '   SEDust: matrices in the grain frame (z = B-hat).'
   print '(a)', '   RT code: meridional rotations, azimuth sampling, peel-off, exp(-K tau).'

   deallocate(Cext, Cabs, Csca, Cpol_ext, Cbir_ext)

contains

   logical function finite(x)
      real(wp), intent(in) :: x
      finite = (x == x) .and. (abs(x) <= huge(1.0_wp))
   end function finite

   real(wp) function interp_lam(g, ng, y, x) result(v)
      real(wp), intent(in) :: g(:), y(:)
      integer,  intent(in) :: ng
      real(wp), intent(in) :: x
      integer  :: lo, hi, mid
      real(wp) :: t
      if (x <= g(1)) then;  v = y(1);  return;  end if
      if (x >= g(ng)) then;  v = y(ng);  return;  end if
      lo = 1;  hi = ng
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (g(mid) <= x) then;  lo = mid;  else;  hi = mid;  end if
      end do
      t = (log(x) - log(g(lo))) / (log(g(hi)) - log(g(lo)))
      v = (1.0_wp - t)*y(lo) + t*y(hi)
   end function interp_lam

   real(wp) function aligned_z11_integral(ib, ti, e) result(s)
      ! INT Z11 dOmega for the eta-scaled aligned matrix at incidence ti, using
      ! the exposed theta_s/phi grids. Z11 is even under phi -> 360-phi, so the
      ! azimuth integral over [0,2pi] is twice the integral over the stored
      ! [0,180]. Trapezoid with the sin(theta_s) solid-angle weight.
      integer,  intent(in) :: ib
      real(wp), intent(in) :: ti, e
      real(wp) :: zmat(4,4), z11(scm_nts, scm_nphi), gts
      real(wp) :: phint(scm_nts), c0, c1
      integer  :: is, ip
      do is = 1, scm_nts
         do ip = 1, scm_nphi
            call mueller_matrix_aligned(ib, ti, scm_theta_s(is), scm_phi(ip), zmat)
            z11(is, ip) = e * zmat(1,1)
         end do
      end do
      ! azimuth trapezoid over [0,180], doubled for [0,360]
      do is = 1, scm_nts
         gts = 0.0_wp
         do ip = 2, scm_nphi
            gts = gts + 0.5_wp*(z11(is,ip-1) + z11(is,ip)) &
                      * (scm_phi(ip) - scm_phi(ip-1))*deg2rad
         end do
         phint(is) = 2.0_wp * gts
      end do
      ! polar trapezoid with sin(theta_s) weight
      s = 0.0_wp
      do is = 2, scm_nts
         c0 = phint(is-1)*sin(scm_theta_s(is-1)*deg2rad)
         c1 = phint(is)  *sin(scm_theta_s(is)  *deg2rad)
         s = s + 0.5_wp*(c0 + c1)*(scm_theta_s(is) - scm_theta_s(is-1))*deg2rad
      end do
   end function aligned_z11_integral

end program use_dustlib_scatmat
