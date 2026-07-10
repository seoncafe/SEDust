program calc_kext_dl07
   !====================================================================
   ! Size-distribution-integrated extinction cross section per H,
   ! albedo, and scattering asymmetry <cos> for the Draine & Li (2007)
   ! / WD01 Milky Way dust model, computed from the SAME optical model
   ! used by build_dl07, and compared point-by-point against Draine's
   ! published table
   !    opt/extcurvs/kext_albedo_WD_MW_3.1_60_D03.all .
   !
   ! Model components (identical to sed_init_dl07):
   !   silicate     : Mie from D03 astrosilicate (q_silicate_full)
   !                  -> Qext, Qsca, Qabs, g.
   !   carbonaceous : ABSORPTION from the DL07 PAH<->graphite xi-blend
   !                  (qpah_dl07, neutral/cation mixed by pah_ionfrac);
   !                  SCATTERING and g from random-oriented D03 graphite
   !                  (q_graphite_full, 1/3 || + 2/3 _|_).  PAH scattering
   !                  is negligible (Rayleigh), so graphite carries it --
   !                  the standard Draine/DL07 treatment.
   !
   ! Per-H integrals over the WD01 size distribution (index 7 = MW R_V=3.1,
   ! b_C=6e-5, the "60" model), with the Draine 2003 x0.93 abundance
   ! reduction ON to match the reference file:
   !   C_X/H   = sum_a (dn/da) pi a^2 Q_X  a dln(a)        [cm^2 / H]
   !   albedo  = C_sca / C_ext
   !   <cos>   = (sum_a (dn/da) pi a^2 Q_sca g a dlna) / C_sca
   !
   ! Output: output/kext_albedo_dl07_ours.dat with our columns alongside
   ! the reference, ready for the comparison plot.
   !====================================================================
   use constants,      only: wp, pi
   use q_silicate_mod, only: q_silicate_full
   use q_graphite_mod, only: q_graphite_full
   use qpah,           only: qpah_dl07, nc_coeff, nc_integer, qpah_use_d03_graphite
   use pah_ioniz_mod,  only: pah_ionfrac
   use grain_dist_mod, only: grain_dist_dl07, gd_apply_d03_reduction
   implicit none

   ! Reference = Draine's ORIGINAL 2003 (Dec 17 2003) kext_albedo table -- the
   ! one his IDL dust_cross() reads and that Qcross_DL07.txt matched.  NB: a
   ! LATER 2009-10-04 recomputation (the plain ".all" in opt/extcurvs/) revised
   ! the FIR opacity up ~12% and the 2175A bump down ~16%; our D03/DL01 optics
   ! match the 2003 table, not the 2009 one.
   character(len=*), parameter :: F_REF = &
      '../data/release/kext_albedo_WD_MW_3.1_60_D03.all_2003'
   character(len=*), parameter :: F_OUT = 'output/kext_albedo_dl07_ours.dat'
   integer,  parameter :: SD_INDEX = 7          ! MW R_V=3.1, b_C=6e-5
   real(wp), parameter :: U_ISRF   = 1.0_wp     ! MMP83 diffuse ISM
   real(wp), parameter :: UM2CM    = 1.0e-4_wp
   real(wp), parameter :: A100     = 0.99999e-2_wp   ! 100 A charge cutoff [um]
   real(wp), parameter :: LAM_MIN  = 1.0e-2_wp  ! restrict to >=100 A (skip X-ray)
   real(wp), parameter :: MDUST_H  = 1.870e-26_wp    ! g dust/H (2003 ref header)

   ! size grid
   integer,  parameter :: NA = 400
   real(wp), parameter :: AMIN = 3.5e-4_wp, AMAX = 3.0_wp   ! um

   real(wp) :: a_um(NA), acm(NA), lna(NA), twt(NA)
   real(wp) :: dn_sil(NA), dn_car(NA), fion(NA)
   real(wp) :: dlna
   integer  :: ja, u, ios, nref, k

   ! reference table
   integer,  parameter :: MAXREF = 2000
   real(wp) :: lam_r(MAXREF), alb_r(MAXREF), g_r(MAXREF), cext_r(MAXREF), kabs_r(MAXREF)
   real(wp) :: dummy

   real(wp) :: lam, Cext, Csca, Cabs, gnum
   real(wp) :: Qext_s, Qsca_s, Qabs_s, g_s
   real(wp) :: Qext_g, Qsca_g, Qabs_g, g_g
   real(wp) :: Qabs_neu, Qabs_ion, Qabs_c, Qsca_c, g_c
   real(wp) :: w

   ! ---- configure the model exactly as build_dl07 does -----------------
   nc_coeff             = 470.0_wp
   nc_integer           = .true.
   qpah_use_d03_graphite = .true.
   gd_apply_d03_reduction = .true.     ! match the reference's x0.93 reduction

   ! ---- size grid (log) + trapezoid-in-log weights ---------------------
   dlna = (log(AMAX) - log(AMIN)) / real(NA-1, wp)
   do ja = 1, NA
      lna(ja) = log(AMIN) + dlna*real(ja-1, wp)
      a_um(ja) = exp(lna(ja))
      acm(ja)  = a_um(ja) * UM2CM
      twt(ja)  = dlna
   end do
   twt(1)  = 0.5_wp*dlna
   twt(NA) = 0.5_wp*dlna

   ! ---- size-distribution weights (wavelength independent) -------------
   do ja = 1, NA
      dn_sil(ja) = grain_dist_dl07(SD_INDEX, 'sil', a_um(ja))
      dn_car(ja) = grain_dist_dl07(SD_INDEX, 'gra', a_um(ja)) &
                 + grain_dist_dl07(SD_INDEX, 'pah', a_um(ja))
      fion(ja)   = pah_ionfrac(a_um(ja), U_ISRF)
      if (a_um(ja) > A100) fion(ja) = 1.0_wp     ! Draine's 100 A cutoff
   end do

   ! ---- read Draine's reference table (lambda grid + ref columns) ------
   open(newunit=u, file=F_REF, status='old', action='read')
   do k = 1, 80
      read(u, '(a)')
   end do
   nref = 0
   do
      ! cols: lambda  albedo  <cos>  C_ext/H  K_abs  <cos^2> [+ "NNNN eV"]
      ! list-directed read takes the first 6 reals, ignores any trailing tokens.
      read(u, *, iostat=ios) lam, alb_r(nref+1), g_r(nref+1), &
                             cext_r(nref+1), kabs_r(nref+1), dummy
      if (ios /= 0) exit
      nref = nref + 1
      lam_r(nref) = lam
      if (nref >= MAXREF) stop 'calc_kext_dl07: MAXREF too small'
   end do
   close(u)
   write(*,'(a,i0,a)') ' read ', nref, ' reference rows'

   ! ---- compute our model at each reference wavelength -----------------
   open(newunit=u, file=F_OUT, status='replace', action='write')
   write(u,'(a)') '# DL07 (WD01 MW R_V=3.1 b_C=6e-5) cross sections vs Draine D03'
   write(u,'(a)') '# computed from build_dl07 optics: sil Mie + carb(qpah abs + graphite sca)'
   write(u,'(a)') '# lambda[um]  Cext_ours[cm2/H]  alb_ours  g_ours  Cabs_ours[cm2/H]' // &
                  '  Cext_ref  alb_ref  g_ref  Cabs_ref[cm2/H]'
   do k = 1, nref
      lam = lam_r(k)
      if (lam < LAM_MIN) cycle
      Cext = 0.0_wp;  Csca = 0.0_wp;  Cabs = 0.0_wp;  gnum = 0.0_wp
      do ja = 1, NA
         ! silicate
         call q_silicate_full(a_um(ja), lam, Qext_s, Qsca_s, Qabs_s, g_s)
         ! carbonaceous: absorption = DL07 xi-blend (charge-mixed);
         !               scattering + g = random-oriented graphite
         call q_graphite_full(a_um(ja), lam, Qext_g, Qsca_g, Qabs_g, g_g)
         call qpah_dl07(0, a_um(ja), lam, Qabs_neu)
         call qpah_dl07(1, a_um(ja), lam, Qabs_ion)
         Qabs_c = (1.0_wp - fion(ja))*Qabs_neu + fion(ja)*Qabs_ion
         Qsca_c = Qsca_g
         g_c    = g_g

         w = pi * acm(ja)*acm(ja) * acm(ja) * twt(ja)   ! pi a^2 * a dlna [cm^3]
         ! extinction = absorption + scattering, per component
         Cext = Cext + w*( dn_sil(ja)*(Qabs_s + Qsca_s) &
                         +  dn_car(ja)*(Qabs_c + Qsca_c) )
         Csca = Csca + w*( dn_sil(ja)*Qsca_s + dn_car(ja)*Qsca_c )
         Cabs = Cabs + w*( dn_sil(ja)*Qabs_s + dn_car(ja)*Qabs_c )
         gnum = gnum + w*( dn_sil(ja)*Qsca_s*g_s + dn_car(ja)*Qsca_c*g_c )
      end do

      block
         real(wp) :: alb, gbar, cabs_ref
         alb  = 0.0_wp;  gbar = 0.0_wp
         if (Cext > 0.0_wp) alb  = Csca / Cext
         if (Csca > 0.0_wp) gbar = gnum / Csca
         cabs_ref = kabs_r(k) * MDUST_H
         write(u,'(es12.5,4(1x,es13.6),1x,es13.6,2(1x,f8.5),1x,es13.6)') &
            lam, Cext, alb, gbar, Cabs, cext_r(k), alb_r(k), g_r(k), cabs_ref
      end block
   end do
   close(u)
   write(*,'(a,a)') ' wrote ', F_OUT
   write(*,'(a)')   ' calc_kext_dl07: done.'

end program calc_kext_dl07
