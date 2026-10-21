program test_absorption_balance
   ! Energy closure of calc_P for one grain.  The GD89 recursion makes the
   ! steady state emit, through its downward rates, exactly the power its
   ! upward rates absorb:
   !    4 pi sum_(u>=2) P(u) kappB(u) = sum_(l<nT) P(l) sum_(u>l) A(u,l) (H_u - H_l),
   ! with the intrabin absorptions of Draine & Li (2001) eq. (28) inside
   ! A(l+1,l).  With the bin-overlap rates every row l < nT absorbs the
   ! support integral 4 pi radiation_absorption, so the emitted power must
   ! equal (1 - P(nT)) times it to rounding, whatever the bin widths are
   ! relative to the photon energies.  The cases take bins from 5e-10 of the
   ! softest field photon to 700 times the hardest, three
   ! field shapes (a narrow band; a broad band with interior zero samples; a
   ! band whose exterior zeros must not extend it), two enthalpy resolutions,
   ! and the wavelength grid both matching and not matching the p_sub_setup
   ! cache.
   use constants, only: wp
   use p_sub, only: p_sub_setup, calc_P, radiation_absorption
   implicit none
   integer, parameter :: NLAM = 400, NFIELD = 3, NGRAIN = 6
   real(wp), parameter :: PI = 3.141592653589793238462643383279502884197_wp
   real(wp), parameter :: HC = 6.62606957e-27_wp*2.99792458e10_wp*1e4_wp   ! erg um
   real(wp), parameter :: TOL = 1e-14_wp
   ! H(T) = H_SCALE T^3 / (1 + (T/300 K)^2) [erg], kappB = K T^6, on a
   ! logarithmic temperature grid [T_LO, T_HI].
   real(wp), parameter :: H_SCALE(NGRAIN) = [3e-20_wp, 3e-18_wp, 3e-17_wp, 3e-16_wp, 3e-15_wp, 3e-14_wp]
   real(wp), parameter :: T_LO(NGRAIN) = [1.0_wp, 3.0_wp, 3.0_wp, 5.0_wp, 10.0_wp, 19.0_wp]
   real(wp), parameter :: T_HI(NGRAIN) = [3000.0_wp, 3000.0_wp, 2000.0_wp, 1000.0_wp, 300.0_wp, 21.0_wp]
   integer, parameter :: NT_CASES(2) = [40, 200]
   real(wp) :: lam(NLAM), lam_use(NLAM), cabs(NLAM), field(NLAM, NFIELD), j_use(NLAM)
   real(wp) :: energy_ev, absorbed, emitted, rel, worst, ratio_lo, ratio_hi, e_lo, e_hi
   real(wp), allocatable :: t(:), h(:), kappb(:), p(:), lnp(:)
   integer :: k, i, f, n, nt, g, ncase, last_k

   do k = 1, NLAM
      lam(k) = 0.05_wp * (2.0e4_wp)**(real(k-1, wp)/real(NLAM-1, wp))
   end do
   cabs = 1e-21_wp * (0.1_wp/lam) * (1.0_wp + 0.5_wp*exp(-((lam - 0.2175_wp)/0.02_wp)**2))
   field = 0.0_wp
   do k = 1, NLAM
      energy_ev = 1.23984_wp / lam(k)
      if (energy_ev >= 4.78_wp .and. energy_ev < 13.6_wp) field(k, 1) = lam(k)**(-2)
      if (energy_ev >= 0.1_wp .and. energy_ev <= 15.0_wp) &
         field(k, 2) = lam(k)**(-1) * exp(-0.3_wp*lam(k))
      if (energy_ev > 2.8_wp .and. energy_ev < 3.2_wp) field(k, 2) = 0.0_wp
      if (energy_ev >= 2.0_wp .and. energy_ev <= 6.0_wp) field(k, 3) = 1.0_wp
   end do
   call p_sub_setup(lam)

   worst = 0.0_wp
   ratio_lo = huge(1.0_wp)
   ratio_hi = 0.0_wp
   ncase = 0
   write(*,'(a)') '  NT grain field grid   min dH/E_soft   max dH/E_hard        P(nT)  emitted/[(1-P(nT)) absorbed]-1'
   do n = 1, size(NT_CASES)
      nt = NT_CASES(n)
      allocate(t(nt), h(nt), kappb(nt), p(nt), lnp(nt))
      do i = 1, NGRAIN
         do k = 1, nt
            t(k) = T_LO(i) * (T_HI(i)/T_LO(i))**(real(k-1, wp)/real(nt-1, wp))
         end do
         h = H_SCALE(i) * t**3 / (1.0_wp + (t/300.0_wp)**2)
         do f = 1, NFIELD
            do g = 1, 2
               lam_use = lam
               if (g == 2) lam_use = 1.001_wp * lam     ! not the cached grid
               ! absorbed power equal to the cooling at the grid's lower third
               kappb = t**6
               j_use = field(:, f) * kappb(nt/3) / radiation_absorption(lam_use, cabs, field(:, f))
               call calc_P(lam_use, cabs, j_use, t, kappb, h, p, lnp)
               absorbed = 4.0_wp * PI * radiation_absorption(lam_use, cabs, j_use)
               emitted  = 4.0_wp * PI * sum(p(2:nt) * kappb(2:nt))
               rel = emitted / ((1.0_wp - p(nt)) * absorbed) - 1.0_wp
               worst = max(worst, abs(rel))
               ! bin widths against the field's photon energies
               last_k = findloc(j_use > 0.0_wp, .true., dim=1, back=.true.)
               e_lo = HC / lam_use(last_k)
               e_hi = HC / lam_use(findloc(j_use > 0.0_wp, .true., dim=1))
               ratio_lo = min(ratio_lo, minval(h(2:nt) - h(1:nt-1)) / e_lo)
               ratio_hi = max(ratio_hi, maxval(h(2:nt) - h(1:nt-1)) / e_hi)
               write(*,'(i4,i6,i6,i5,2es16.3,es13.3,es16.3)') nt, i, f, g, &
                  minval(h(2:nt) - h(1:nt-1)) / e_lo, maxval(h(2:nt) - h(1:nt-1)) / e_hi, &
                  p(nt), rel
               ncase = ncase + 1
            end do
         end do
      end do
      deallocate(t, h, kappb, p, lnp)
   end do
   write(*,'(a,i0,a,es10.3,a,es10.3,a,es10.3)') 'cases ', ncase, &
      ', bin width / photon energy from ', ratio_lo, ' to ', ratio_hi, &
      ', worst |closure - 1| ', worst
   if (ratio_lo > 1e-9_wp .or. ratio_hi < 10.0_wp) error stop 'cases do not span the bin widths'
   if (worst > TOL) error stop 'calc_P energy closure'
   write(*,'(a)') 'PASS: calc_P emits the support absorption to rounding'
end program test_absorption_balance
