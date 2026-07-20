program compare_birefringence
   ! Certification of the birefringence (circular-polarization) optic added to
   ! the orientation-resolved pipeline.
   !
   ! The forward-scattering amplitude S_fwd(jori) is computed by AMPL from the
   ! converged T-matrix.  Its IMAGINARY part in the forward direction gives the
   ! extinction, Q_ext(jori) = (4 pi/k) Im[S_fwd]/(pi a_eff^2), and the
   ! polarized extinction is C_pol,ext = 0.5*(Q_ext(3)-Q_ext(2))*pi a_eff^2.
   ! Its REAL part is the birefringence twin,
   !   Q_re(jori) = (4 pi/k) Re[S_fwd]/(pi a_eff^2),
   !   C_bir = 0.5*(Q_re(3)-Q_re(2))*pi a_eff^2,
   ! the phase retardation that converts Stokes U into V on propagation.  This
   ! is the U<->V coupling coefficient (Mishchenko's extinction-matrix element
   ! K_34; Martin 1974 MNRAS 167, 57; Whitney & Wolff 2002 ApJ 574, 205), and
   ! we adopt exactly the orientation and sign bookkeeping already used for
   ! C_pol,ext (jori=3 minus jori=2, halved).
   !
   ! There is no external orientation-resolved birefringence benchmark to
   ! compare against, so the certification is an INTERNAL Kramers-Kronig test.
   ! The forward-amplitude difference between the two linear polarizations,
   !   D(omega) = S_fwd(jori=3) - S_fwd(jori=2),
   ! is analytic in the upper-half omega-plane (causality) and vanishes as
   ! omega -> infinity (the extinction paradox is polarization-independent, so
   ! Q_ext(2) = Q_ext(3) at geometric-optics order, i.e. Im[D] and Re[D] -> 0).
   ! With the crossing relation D(-omega) = D*(omega) the unsubtracted
   ! one-sided Kramers-Kronig relation therefore holds directly for D:
   !
   !   Re[D](omega) = (2/pi) P INT_0^inf omega' Im[D](omega')/(omega'^2-omega^2) domega'.
   !
   ! CAUTION applied here (the omega-power that made this work): the KK pair is
   ! the (m-1)-like difference dm(omega), NOT Q and NOT the bare amplitude D.
   ! With the code's optical theorem C_ext = (4 pi/k) Im[S], the stored S is the
   ! amplitude f, so the effective-medium refractive index is m-1 ~ S/k^2 and
   !   dm = m3 - m2 ~ (C_bir + i C_pol) / omega        (omega = 1/lambda),
   ! since the cross section is C = Q * area ~ omega * dm.  In omega-weight:
   ! dm (KK-clean, omega^0) -> C, Q (one extra omega) -> D = S3-S2 (two extra).
   ! Feeding Q (extra cext_fac = 2 lambda) or the bare amplitude D to the naive
   ! relation therefore fails; feeding dm = C/omega succeeds.  The omega power
   ! is pinned analytically by the Rayleigh limit, where the static
   ! polarizability makes Im[dm] -> const while C_pol ~ omega as omega -> 0.
   ! C_pol = 0.5*(Q_ext(3)-Q_ext(2))*pi a_eff^2 is the validated dichroism (few
   ! parts in 1e-4 vs HD23) and C_bir = 0.5*(Q_re(3)-Q_re(2))*pi a_eff^2 the
   ! new birefringence, both from the same tmatrix_oriented_cross outputs (no
   ! separate AMPL call); agreement of C_bir with the Hilbert transform of
   ! C_pol certifies the real-part optic regardless of the overall sign choice.
   !
   ! Because the sampled band is finite, the KK tails are truncated and the
   ! edges are spoiled; the residual is reported over the interior (middle 80%)
   ! of the band, per the standard practice for finite-range Hilbert pairs.
   !
   ! The program also runs:
   !   - the Rayleigh-analytic anchor: the closed-form small-x birefringence
   !     0.5*(8 pi/(lam a^2))*(Re alpha_b - Re alpha_a) vs the pipeline value;
   !   - the continuity anchor: T-matrix birefringence just above x=0.1 vs the
   !     Rayleigh value just below, at fixed radius;
   !   - the precision anchor: does C_bir recovered from the 5/6-significant-
   !     figure stored Q_re match the directly computed C_bir, or does the
   !     Q_re(3)-Q_re(2) cancellation demand storing the difference directly.

   use, intrinsic :: iso_fortran_env, only: error_unit
   use constants,        only: wp
   use read_index,       only: load_index, interp_m
   use tmatrix_oriented, only: tmatrix_oriented_cross, oriented_cross_sections
   use asymptotic_optics, only: rayleigh_limit
   implicit none

   character(len=*), parameter :: f_aeff  = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: f_wave  = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: f_index = &
      '../data/dielectric/index_DH21Ad_P0.20_0.00_1.400'

   integer,  parameter :: NA = 169, NW = 1129
   real(wp), parameter :: EPS_BA  = 1.4_wp
   real(wp), parameter :: DDELT   = 1.0e-3_wp
   integer,  parameter :: NDGS    = 2, NP_OBL = -1
   real(wp), parameter :: X_SMALL = 0.1_wp, X_LARGE = 50.0_wp
   real(wp), parameter :: PI      = acos(-1.0_wp)

   real(wp) :: a_eff(NA), lambda(NW)
   complex(wp) :: m_cache(NW)
   real(wp) :: nr, ki
   integer  :: jw

   call read_one_col(f_aeff, NA, a_eff)
   call read_one_col(f_wave, NW, lambda)
   call load_index(f_index)
   do jw = 1, NW
      call interp_m(lambda(jw), nr, ki)
      m_cache(jw) = cmplx(nr, ki, kind=wp)
   end do

   write(*,'(a)') '======================================================================'
   write(*,'(a)') ' Birefringence certification (real part of forward-amplitude diff)'
   write(*,'(a)') '======================================================================'

   call kk_check()
   call rayleigh_analytic_check()
   call continuity_check()
   call precision_check()

contains

   ! ------------------------------------------------------------------
   ! Kramers-Kronig certification: Re[D] vs Hilbert transform of Im[D].
   ! ------------------------------------------------------------------
   subroutine kk_check()
      integer, parameter :: NR_TEST = 4
      integer  :: rad_idx(NR_TEST) = (/ 60, 80, 100, 120 /)   ! grid radii
      integer  :: ir, ia, stride, n, i
      real(wp), allocatable :: om(:), red(:), imd(:), red_kk(:), reg(:), img(:)
      real(wp) :: med_rel, max_scale

      ! First validate the discrete transform itself on a Lorentz oscillator,
      ! whose Re/Im are an exact KK pair, so any residual below is physics or
      ! finite-band truncation, not the transform.
      call kk_selftest()

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Kramers-Kronig: C_bir from Hilbert transform of C_pol'
      write(*,'(a)') '   KK pair is the (m-1)-like difference dm(om) ~ (C_bir + i C_pol)*lam'
      write(*,'(a)') '   = (C_bir + i C_pol)/om.  This is NOT Q (extra cext_fac=2*lam) and'
      write(*,'(a)') '   NOT the amplitude D=S3-S2 (extra om^2): with C_ext=(4pi/k)Im[S] the'
      write(*,'(a)') '   code S is the amplitude f, so m-1 ~ S/k^2 and dm ~ C/om.  C_pol='
      write(*,'(a)') '   0.5(Qext3-Qext2)*area (dichroism, imag), C_bir=0.5(Qre3-Qre2)*area'
      write(*,'(a)') '   (birefringence, real); the /om power is pinned by the Rayleigh limit'
      write(*,'(a)') '   (Im[dm]->const while C_pol~om as om->0).  om=1/lam, subtractive KK,'
      write(*,'(a)') '   FULL grid (Rayleigh+T-matrix), interior middle 80%.'
      write(*,'(a)') '   The grid stops at lambda=0.0912 um (912 A), where the dichroism is'
      write(*,'(a)') '   still strong, so the KK high-om tail is truncated; a missing distant'
      write(*,'(a)') '   band adds a near-constant offset to Re, removed by once-subtracting'
      write(*,'(a)') '   (anchoring the reconstruction at the interior peak of C_bir).'
      write(*,'(a)') '   a_eff[um]  n_grid  median|Cbir_kk-Cbir|/max   median|Cbir_kk/Cbir-1|'

      ! The KK integral needs the full omega axis: C_pol is largest in the
      ! T-matrix band and vanishes toward both ends.  Gather the pair across the
      ! entire wavelength grid through the regime dispatcher.  The transform is
      ! evaluated directly on the native (log-spaced, hence nonuniform) om grid
      ! by singularity subtraction, so no resampling artifact is introduced.
      stride = max(1, NW/1100)               ! keep essentially the full grid
      allocate(om(NW), red(NW), imd(NW), red_kk(NW), reg(NW), img(NW))
      do ir = 1, NR_TEST
         ia = rad_idx(ir)
         call gather_amplitude(ia, stride, n, om, red, imd)
         if (n < 20) then
            write(*,'(f10.5,i9,a)') a_eff(ia), n, '   grid too small, skipped'
            cycle
         end if
         call sort3(n, om, red, imd)
         ! KK-clean (m-1)-like pair is dm ~ (C_bir + i C_pol)*lambda, i.e. the
         ! cross section divided by om (verified analytically in the Rayleigh
         ! limit: Im[dm] -> const while C_pol ~ om as om -> 0).
         do i = 1, n
            reg(i) = red(i) / om(i)
            img(i) = imd(i) / om(i)
         end do
         call subtractive_kk(n, om, img, red_kk)
         call anchor_once_subtracted(n, reg, red_kk)
         call interior_residual(n, reg, red_kk, med_rel, max_scale)
         write(*,'(f10.5,i9,2(es22.4))') a_eff(ia), n, &
              interior_median_abs(n, reg, red_kk)/max_scale, med_rel
      end do
      deallocate(om, red, imd, red_kk, reg, img)
   end subroutine kk_check


   subroutine kk_selftest()
      ! Validate the singularity-subtracted one-sided KK on a Lorentz oscillator
      !   chi(om) = om_p^2 / (om0^2 - om^2 - i*gamma*om),
      ! Re even, Im odd, chi -> 0 at infinity: an exact one-sided KK pair.  Run
      ! on a log-spaced om grid to mimic the native (1/lambda) sampling.
      integer, parameter :: NU = 1200
      real(wp) :: u(NU), imc(NU), rec(NU), rekk(NU)
      real(wp), parameter :: om0 = 3.0_wp, gam = 0.6_wp, omp = 1.0_wp
      real(wp), parameter :: umin = 1.0e-3_wp, umax = 3.0e2_wp
      integer  :: i
      real(wp) :: d2, med_rel, max_scale, r
      r = (umax/umin)**(1.0_wp/real(NU-1, wp))
      do i = 1, NU
         u(i)   = umin * r**(i-1)             ! log-spaced, like 1/lambda
         d2     = (om0*om0 - u(i)*u(i))
         imc(i) = omp*omp * gam*u(i) / (d2*d2 + gam*gam*u(i)*u(i))
         rec(i) = omp*omp * d2       / (d2*d2 + gam*gam*u(i)*u(i))
      end do
      call subtractive_kk(NU, u, imc, rekk)
      call interior_residual(NU, rec, rekk, med_rel, max_scale)
      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Transform self-test: subtractive KK on a Lorentz oscillator (log grid)'
      write(*,'(a,es12.4,a,es12.4)') '   interior median |Re_kk/Re-1| = ', med_rel, &
           '   median|Re_kk-Re|/max|Re| = ', interior_median_abs(NU, rec, rekk)/max_scale
   end subroutine kk_selftest


   subroutine anchor_once_subtracted(n, f, g)
      ! Convert the (truncated) full KK reconstruction g into the once-subtracted
      ! result anchored at the interior node where |f| is largest: g <- g - [g-f]
      ! there.  This removes the near-constant Re offset from the missing
      ! dichroism beyond the grid's high-om edge, leaving the shape to be tested.
      integer,  intent(in)    :: n
      real(wp), intent(in)    :: f(:)
      real(wp), intent(inout) :: g(:)
      integer  :: lo, hi, i, ianch
      real(wp) :: fmax, offset
      lo = max(1, int(0.1_wp*n) + 1);  hi = min(n, int(0.9_wp*n))
      fmax = -1.0_wp;  ianch = lo
      do i = lo, hi
         if (abs(f(i)) > fmax) then
            fmax = abs(f(i));  ianch = i
         end if
      end do
      offset = g(ianch) - f(ianch)
      g(1:n) = g(1:n) - offset
   end subroutine anchor_once_subtracted


   subroutine subtractive_kk(n, om, imf, ref)
      ! One-sided Kramers-Kronig by singularity subtraction, evaluated on the
      ! (possibly nonuniform, sorted-ascending) grid om.  Using the identity
      ! P INT_0^inf dom'/(om'^2-om^2) = 0,
      !   Re(om_i) = (2/pi) INT_0^inf [om' Im(om') - om_i Im(om_i)]
      !                                 /(om'^2 - om_i^2) dom',
      ! whose integrand is regular at om' = om_i (removable pole).  The integral
      ! is done by the trapezoid rule; at the removable node the analytic limit
      ! [Im(om_i) + om_i Im'(om_i)]/(2 om_i) is used, with Im' a grid-centered
      ! finite difference.  This needs no uniform resampling and no principal-
      ! value cancellation.
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: om(:), imf(:)
      real(wp), intent(out) :: ref(:)
      integer  :: i, j
      real(wp) :: integ, deriv
      real(wp), allocatable :: h(:)
      allocate(h(n))
      do i = 1, n
         do j = 1, n
            if (j == i) then
               if (i == 1) then
                  deriv = (imf(2) - imf(1)) / (om(2) - om(1))
               else if (i == n) then
                  deriv = (imf(n) - imf(n-1)) / (om(n) - om(n-1))
               else
                  deriv = (imf(i+1) - imf(i-1)) / (om(i+1) - om(i-1))
               end if
               h(j) = (imf(i) + om(i)*deriv) / (2.0_wp*om(i))
            else
               h(j) = (om(j)*imf(j) - om(i)*imf(i)) / (om(j)*om(j) - om(i)*om(i))
            end if
         end do
         integ = 0.0_wp
         do j = 1, n-1
            integ = integ + 0.5_wp*(h(j) + h(j+1))*(om(j+1) - om(j))
         end do
         ref(i) = (2.0_wp/PI) * integ
      end do
      deallocate(h)
   end subroutine subtractive_kk


   real(wp) function interior_median_abs(n, f, g) result(med)
      ! Median |f-g| over the interior middle 80%.
      integer,  intent(in) :: n
      real(wp), intent(in) :: f(:), g(:)
      integer :: lo, hi, m, i
      real(wp), allocatable :: d(:)
      lo = max(1, int(0.1_wp*n) + 1);  hi = min(n, int(0.9_wp*n))
      m  = hi - lo + 1
      allocate(d(m))
      do i = lo, hi
         d(i-lo+1) = abs(f(i) - g(i))
      end do
      call sort1(m, d)
      med = d((m+1)/2)
      deallocate(d)
   end function interior_median_abs


   subroutine interior_residual(n, f, g, med_rel, max_scale)
      ! Interior middle-80% residual statistics comparing g (reconstructed) to
      ! f (direct): max|f| as a robust scale, and the median pointwise relative
      ! deviation where the signal exceeds 0.05*max|f|.
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: f(:), g(:)
      real(wp), intent(out) :: med_rel, max_scale
      integer  :: lo, hi, i, m
      real(wp), allocatable :: dr(:)
      lo = max(1, int(0.1_wp*n) + 1);  hi = min(n, int(0.9_wp*n))
      max_scale = 0.0_wp
      do i = lo, hi
         if (abs(f(i)) > max_scale) max_scale = abs(f(i))
      end do
      if (max_scale <= 0.0_wp) max_scale = 1.0_wp
      allocate(dr(hi-lo+1))
      m = 0
      do i = lo, hi
         if (abs(f(i)) > 0.05_wp*max_scale) then
            m = m + 1
            dr(m) = abs(g(i)/f(i) - 1.0_wp)
         end if
      end do
      if (m > 0) then
         call sort1(m, dr)
         med_rel = dr((m+1)/2)
      else
         med_rel = 0.0_wp
      end if
      deallocate(dr)
   end subroutine interior_residual


   subroutine gather_amplitude(ia, stride, n, om, red, imd)
      ! Fill om (=1/lambda, ~frequency) and the KK pair (C_bir, C_pol) over the
      ! FULL wavelength grid for radius index ia, subsampled by stride.
      ! oriented_cross_sections dispatches Rayleigh / T-matrix / geometric
      ! optics, so the pair is defined at every node and the KK integral spans
      ! the whole omega axis.
      !
      ! The (m-1)-like difference between the two linear polarizations,
      !   dm(om) ~ [S_fwd(3) - S_fwd(2)] / k  ~  C_pol + i C_bir,
      ! is the analytic quantity whose real and imaginary parts are the KK
      ! pair; forming the cross section C = Q * area removes the omega factor
      ! that the amplitude S and the efficiency Q each carry, so no explicit
      ! multiply by lambda is needed (it cancels cext_fac = 2*lambda).
      !   red = C_bir = 0.5*(Qre3-Qre2)*area   (birefringence, real part)
      !   imd = C_pol = 0.5*(Qext3-Qext2)*area (dichroism,     imag part)
      integer,  intent(in)  :: ia, stride
      integer,  intent(out) :: n
      real(wp), intent(out) :: om(:), red(:), imd(:)
      integer  :: j, flag
      real(wp) :: area
      real(wp) :: qext_ori(3), qabs_ori(3), qsca_ori(3), qre_ori(3)

      n = 0
      area = PI * a_eff(ia) * a_eff(ia)
      do j = 1, NW, stride
         call oriented_cross_sections(a_eff(ia), lambda(j), m_cache(j), EPS_BA, &
                  NP_OBL, DDELT, NDGS, qext_ori, qabs_ori, qsca_ori, flag, &
                  qre_ori=qre_ori)
         n = n + 1
         om(n)  = 1.0_wp / lambda(j)
         red(n) = 0.5_wp * (qre_ori(3)  - qre_ori(2))  * area
         imd(n) = 0.5_wp * (qext_ori(3) - qext_ori(2)) * area
      end do
   end subroutine gather_amplitude


   ! ------------------------------------------------------------------
   ! Rayleigh-analytic anchor.
   ! ------------------------------------------------------------------
   subroutine rayleigh_analytic_check()
      integer  :: ia, jw2, ib
      real(wp) :: x, area, qbir_pipe, qbir_ana, rel
      real(wp) :: qext_o(3), qsca_o(3), walb, asymm, qre_o(3)
      logical  :: found

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Rayleigh analytic anchor: pipeline C_bir vs closed form (x < 0.1)'
      write(*,'(a)') '   a_eff[um]  lambda[um]     x       C_bir_pipe     C_bir_ana    rel.dev'

      ! A few small-x nodes.
      do ib = 1, 3
         ia = 30 + 20*ib
         found = .false.
         do jw2 = NW, 1, -1               ! long wavelengths first (small x)
            x = 2.0_wp*PI*a_eff(ia)/lambda(jw2)
            if (x < 0.08_wp .and. x > 0.02_wp) then
               found = .true.
               exit
            end if
         end do
         if (.not. found) cycle
         area = PI * a_eff(ia) * a_eff(ia)
         call rayleigh_limit(a_eff(ia), lambda(jw2), real(m_cache(jw2),wp), &
                 aimag(m_cache(jw2)), EPS_BA, qext_o(1), qsca_o(1), walb, asymm, &
                 qext_ori=qext_o, qabs_ori=qsca_o, qsca_ori=qsca_o, qre_ori=qre_o)
         qbir_pipe = 0.5_wp*(qre_o(3) - qre_o(2)) * area
         call rayleigh_birefringence_closed(a_eff(ia), lambda(jw2), m_cache(jw2), &
                 EPS_BA, qbir_ana)
         rel = abs(qbir_pipe/qbir_ana - 1.0_wp)
         write(*,'(f10.5,f11.4,es11.3,3es14.5)') a_eff(ia), lambda(jw2), x, &
              qbir_pipe, qbir_ana, rel
      end do
   end subroutine rayleigh_analytic_check


   subroutine rayleigh_birefringence_closed(a_eff1, lam, m, eps_ba, cbir)
      ! Independent closed form of the small-x birefringence, derived here from
      ! the spheroid dipole polarizabilities (Draine 1992), NOT via the
      ! pipeline routine, so agreement checks the wiring and the prefactor.
      !   C_bir = 0.5*(8 pi/(lam a^2))*(Re alpha_b - Re alpha_a) * pi a^2
      real(wp),    intent(in)  :: a_eff1, lam, eps_ba
      complex(wp), intent(in)  :: m
      real(wp),    intent(out) :: cbir
      real(wp)    :: axrat, e2, e, ala, alb, fac, area
      complex(wp) :: eps, alpha_a, alpha_b
      real(wp) :: n_r, k_i
      n_r = real(m, wp);  k_i = abs(aimag(m))
      eps = cmplx(n_r*n_r - k_i*k_i, 2.0_wp*n_r*k_i, kind=wp)
      axrat = 1.0_wp / eps_ba
      e2 = abs(1.0_wp - 1.0_wp/(axrat*axrat))
      e  = sqrt(e2)
      if (axrat < 1.0_wp) then
         ala = (1.0_wp + 1.0_wp/e2) * (1.0_wp - atan(e)/e)
      else if (axrat > 1.0_wp) then
         ala = (1.0_wp/e2 - 1.0_wp) * &
               (log((1.0_wp + e)/(1.0_wp - e))/(2.0_wp*e) - 1.0_wp)
      else
         ala = 1.0_wp/3.0_wp
      end if
      alb = (1.0_wp - ala) / 2.0_wp
      fac     = a_eff1**3 / 3.0_wp
      alpha_a = fac * (eps - 1.0_wp) / ((eps - 1.0_wp)*ala + 1.0_wp)
      alpha_b = fac * (eps - 1.0_wp) / ((eps - 1.0_wp)*alb + 1.0_wp)
      fac  = 8.0_wp * PI / (lam * a_eff1*a_eff1)
      area = PI * a_eff1 * a_eff1
      cbir = 0.5_wp * fac * (real(alpha_b,wp) - real(alpha_a,wp)) * area
   end subroutine rayleigh_birefringence_closed


   ! ------------------------------------------------------------------
   ! Continuity anchor across x = 0.1.
   ! ------------------------------------------------------------------
   subroutine continuity_check()
      integer  :: ia, jw_ray, jw_tm, j, ierr
      real(wp) :: x, area, qbir_ray, qbir_tm, rel
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3), qre_o(3), walb, asymm

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Continuity anchor across x=0.1 (Rayleigh below vs T-matrix above)'
      write(*,'(a)') '   a_eff[um]   x_ray    x_tm    Qbir_ray     Qbir_tm      rel.dev'

      do ia = 70, 110, 20
         area = PI * a_eff(ia) * a_eff(ia)
         ! nearest grid wavelength with x just below / just above 0.1
         jw_ray = 0;  jw_tm = 0
         do j = 1, NW
            x = 2.0_wp*PI*a_eff(ia)/lambda(j)
            if (x < X_SMALL) then
               if (jw_ray == 0) jw_ray = j    ! first (smallest lam step) below
            end if
         end do
         ! choose the pair straddling 0.1 most tightly
         do j = 1, NW-1
            x = 2.0_wp*PI*a_eff(ia)/lambda(j)
            if (x >= X_SMALL .and. 2.0_wp*PI*a_eff(ia)/lambda(j+1) < X_SMALL) then
               jw_tm  = j
               jw_ray = j+1
               exit
            end if
         end do
         if (jw_tm == 0 .or. jw_ray == 0) cycle
         call rayleigh_limit(a_eff(ia), lambda(jw_ray), real(m_cache(jw_ray),wp), &
                 aimag(m_cache(jw_ray)), EPS_BA, qext_o(1), qsca_o(1), walb, asymm, &
                 qext_ori=qext_o, qabs_ori=qabs_o, qsca_ori=qsca_o, qre_ori=qre_o)
         qbir_ray = 0.5_wp*(qre_o(3) - qre_o(2))
         call tmatrix_oriented_cross(a_eff(ia), lambda(jw_tm), m_cache(jw_tm), &
                 EPS_BA, NP_OBL, DDELT, NDGS, qext_o, qsca_o, qabs_o, ierr, qre_ori=qre_o)
         if (ierr /= 0) cycle
         qbir_tm = 0.5_wp*(qre_o(3) - qre_o(2))
         if (qbir_ray /= 0.0_wp) then
            rel = abs(qbir_tm/qbir_ray - 1.0_wp)
         else
            rel = 0.0_wp
         end if
         write(*,'(f10.5,2f9.4,3es13.4)') a_eff(ia), &
              2.0_wp*PI*a_eff(ia)/lambda(jw_ray), 2.0_wp*PI*a_eff(ia)/lambda(jw_tm), &
              qbir_ray, qbir_tm, rel
      end do
   end subroutine continuity_check


   ! ------------------------------------------------------------------
   ! Precision anchor: stored 5/6-sig-fig Q_re vs direct C_bir.
   ! ------------------------------------------------------------------
   subroutine precision_check()
      integer  :: ia, j, ierr, ncnt
      real(wp) :: x, area, qbir_dir, qbir_sto, rel, cancel
      real(wp) :: qext_o(3), qsca_o(3), qabs_o(3), qre_o(3)
      real(wp) :: q2rt, q3rt
      real(wp), allocatable :: relv(:), canv(:)
      real(wp) :: med_rel, max_rel, med_can, max_can

      write(*,'(a)') '----------------------------------------------------------------------'
      write(*,'(a)') ' Precision anchor: C_bir from es13.5-stored Q_re vs direct C_bir'
      write(*,'(a)') '   (cancel = |Q_re(2)|/|Q_re(3)-Q_re(2)|, the loss of significance)'

      allocate(relv(NW*3), canv(NW*3))
      ncnt = 0
      do ia = 40, 130, 10
         area = PI * a_eff(ia) * a_eff(ia)
         do j = 1, NW, 20
            x = 2.0_wp*PI*a_eff(ia)/lambda(j)
            if (x <= X_SMALL .or. x >= X_LARGE) cycle
            call tmatrix_oriented_cross(a_eff(ia), lambda(j), m_cache(j), EPS_BA, &
                     NP_OBL, DDELT, NDGS, qext_o, qsca_o, qabs_o, ierr, qre_ori=qre_o)
            if (ierr /= 0) cycle
            qbir_dir = 0.5_wp*(qre_o(3) - qre_o(2)) * area
            q2rt = roundtrip_es135(qre_o(2))
            q3rt = roundtrip_es135(qre_o(3))
            qbir_sto = 0.5_wp*(q3rt - q2rt) * area
            if (qbir_dir == 0.0_wp) cycle
            ncnt = ncnt + 1
            relv(ncnt) = abs(qbir_sto/qbir_dir - 1.0_wp)
            if (abs(qre_o(3)-qre_o(2)) > 0.0_wp) then
               canv(ncnt) = abs(qre_o(2)) / abs(qre_o(3)-qre_o(2))
            else
               canv(ncnt) = huge(1.0_wp)
            end if
         end do
      end do

      if (ncnt > 0) then
         call sort1(ncnt, relv)
         med_rel = relv((ncnt+1)/2);  max_rel = relv(ncnt)
         call sort1(ncnt, canv)
         med_can = canv((ncnt+1)/2);  max_can = canv(ncnt)
         write(*,'(a,i0)')       '   nodes tested            : ', ncnt
         write(*,'(a,es12.4)')   '   median rel. err (stored): ', med_rel
         write(*,'(a,es12.4)')   '   max    rel. err (stored): ', max_rel
         write(*,'(a,es12.4)')   '   median cancellation     : ', med_can
         write(*,'(a,es12.4)')   '   max    cancellation     : ', max_can
         if (max_rel > 1.0e-2_wp) then
            write(*,'(a)') '   VERDICT: per-orientation Q_re at es13.5 loses the birefringence;'
            write(*,'(a)') '            store the difference 0.5*(Q_re3-Q_re2) directly instead.'
         else
            write(*,'(a)') '   VERDICT: es13.5 per-orientation Q_re recovers C_bir adequately.'
         end if
      end if
      deallocate(relv, canv)
   end subroutine precision_check


   real(wp) function roundtrip_es135(v) result(r)
      ! Round v through the ES13.5 stored representation (the write format used
      ! by run_q_jori.f90) and read it back.
      real(wp), intent(in) :: v
      character(len=16) :: s
      write(s,'(es13.5)') v
      read(s,*) r
   end function roundtrip_es135


   ! ------------------------------------------------------------------
   ! utilities
   ! ------------------------------------------------------------------
   subroutine sort1(n, a)
      integer,  intent(in)    :: n
      real(wp), intent(inout) :: a(:)
      integer :: i, j
      real(wp) :: k
      do i = 2, n
         k = a(i);  j = i-1
         do while (j >= 1)
            if (a(j) <= k) exit
            a(j+1) = a(j);  j = j-1
         end do
         a(j+1) = k
      end do
   end subroutine sort1

   subroutine sort3(n, a, b, c)
      ! Sort a ascending, carrying b and c.
      integer,  intent(in)    :: n
      real(wp), intent(inout) :: a(:), b(:), c(:)
      integer :: i, j
      real(wp) :: ka, kb, kc
      do i = 2, n
         ka = a(i);  kb = b(i);  kc = c(i);  j = i-1
         do while (j >= 1)
            if (a(j) <= ka) exit
            a(j+1) = a(j);  b(j+1) = b(j);  c(j+1) = c(j);  j = j-1
         end do
         a(j+1) = ka;  b(j+1) = kb;  c(j+1) = kc
      end do
   end subroutine sort3

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

end program compare_birefringence
