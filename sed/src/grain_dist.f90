module grain_dist_mod
   !--------------------------------------------------------------------
   ! Weingartner & Draine (2001) / Draine & Li (2007) grain size
   ! distributions for silicate, graphite, and PAH grains.
   !
   ! Ported from Qext_new/grain_dist.f90 (GRAIN_DIST_DL07), originally
   ! from B.T. Draine's web distribution (http://physics.gmu.edu/~joe/),
   ! with the Draine & Li (2007) VSG update and the LMC/SMC parameter
   ! sets added by K.-I. Seon. Converted here to a free-form F90 module
   ! using the project working precision `wp`.
   !
   !   grain_dist_dl07(index, dtype, a_um) -> dn/da per H  [cm^-1]
   !     dtype = 'sil' (silicate), 'gra' (graphite), 'pah' (PAH);
   !             matched case-insensitively on the first 3 characters.
   !     Full carbonaceous dn/da = 'gra' + 'pah' (the xi_gra split below
   !     partitions the carbonaceous distribution into the two).
   !
   !   index -> (R_V, 10^5 b_C, case):
   !     1-7   : R_V=3.1, b_C=0..6   (case A)   [7 = MW3.1 q_PAH=4.58%]
   !     8-12  : R_V=4.0, b_C=0..4   (case A)
   !     13-16 : R_V=5.5, b_C=0..3   (case A)
   !     17-21 : R_V=4.0, b_C=0..4   (case B)
   !     22-25 : R_V=5.5, b_C=0..3   (case B)
   !     26-28 : LMC avg, b_C=0..2
   !     29-31 : LMC 2,   b_C=0,0.5,1.0
   !     32    : SMC bar, b_C=0
   !--------------------------------------------------------------------
   use constants, only: wp
   implicit none
   private
   public :: grain_dist_dl07
   public :: gd_apply_d03_reduction

   ! Whether to apply the Draine (2003, ARA&A) *0.93 abundance reduction to
   ! the MW (index<=25) distributions. Qext_new/grain_dist.f90 applies it,
   ! but the DL07 published spectra
   ! use the un-reduced WD01 abundances. Default .false. to match DL07.
   logical, save :: gd_apply_d03_reduction = .false.

   ! WD01/DL07 fit parameters (index 1..32), from Qext_new/grain_dist.f90
   real(wp), parameter :: ALPHAGARR(32) = [ &
      -2.25_wp,-2.17_wp,-2.04_wp,-1.91_wp,-1.84_wp,-1.72_wp,-1.54_wp,-2.26_wp,-2.16_wp,-2.01_wp, &
      -1.83_wp,-1.64_wp,-2.35_wp,-2.12_wp,-1.94_wp,-1.61_wp,-2.62_wp,-2.52_wp,-2.36_wp,-2.09_wp, &
      -1.96_wp,-2.80_wp,-2.67_wp,-2.45_wp,-1.90_wp, &
      -2.91_wp,-2.99_wp, 4.43_wp,-2.94_wp,-2.82_wp, 4.16_wp,-2.79_wp]
   real(wp), parameter :: BETAGARR(32) = [ &
      -0.0648_wp,-0.0382_wp,-0.111_wp,-0.125_wp,-0.132_wp,-0.322_wp,-0.165_wp,-0.199_wp,-0.0862_wp, &
      -0.0973_wp,-0.175_wp,-0.247_wp,-0.668_wp,-0.67_wp,-0.853_wp,-0.722_wp,-0.0144_wp,-0.0541_wp, &
      -0.0957_wp,-0.193_wp,-0.813_wp,0.0356_wp,0.0129_wp,-0.00132_wp,-0.0517_wp, &
       0.895_wp, 2.460_wp, 0.000_wp, 5.220_wp, 9.010_wp, 0.000_wp, 1.120_wp]
   real(wp), parameter :: ATGARR(32) = [ &
      0.00745_wp,0.00373_wp,0.00828_wp,0.00837_wp,0.00898_wp,0.0254_wp,0.0107_wp,0.0241_wp, &
      0.00867_wp,0.00811_wp,0.0117_wp,0.0152_wp,0.148_wp,0.0686_wp,0.0786_wp,0.0418_wp,0.0187_wp, &
      0.0366_wp,0.0305_wp,0.0199_wp,0.0693_wp,0.0203_wp,0.0134_wp,0.0275_wp,0.012_wp, &
      0.5780_wp,0.0980_wp,0.00322_wp,0.3730_wp,0.3920_wp,0.3420_wp,0.0190_wp]
   real(wp), parameter :: ACGARR(32) = [ &
      0.606_wp,0.586_wp,0.543_wp,0.499_wp,0.489_wp,0.438_wp,0.428_wp,0.861_wp,0.803_wp,0.696_wp, &
      0.604_wp,0.536_wp,1.96_wp,1.35_wp,0.921_wp,0.72_wp,5.74_wp,6.65_wp,6.44_wp,4.6_wp,3.48_wp,3.43_wp, &
      3.44_wp,5.14_wp,7.28_wp, &
      1.210_wp,0.6410_wp,0.2850_wp,0.3490_wp,0.2690_wp,0.04930_wp,0.5220_wp]
   real(wp), parameter :: CGARR(32) = [ &
      9.94e-11_wp,3.79e-10_wp,5.57e-11_wp,4.15e-11_wp,2.90e-11_wp,3.20e-12_wp,9.99e-12_wp, &
      5.47e-12_wp,4.58e-11_wp,3.96e-11_wp,1.42e-11_wp,5.83e-12_wp,4.82e-14_wp,3.65e-13_wp, &
      2.57e-13_wp,7.58e-13_wp,6.46e-12_wp,1.08e-12_wp,1.62e-12_wp,4.21e-12_wp,2.95e-13_wp, &
      2.74e-12_wp,7.25e-12_wp,8.79e-13_wp,2.86e-12_wp, &
      7.12e-17_wp,3.51e-15_wp,9.57e-24_wp,9.92e-17_wp,6.20e-17_wp,3.05e-15_wp,8.36e-14_wp]
   real(wp), parameter :: ALPHASARR(32) = [ &
      -1.48_wp,-1.46_wp,-1.43_wp,-1.41_wp,-2.1_wp,-2.1_wp,-2.21_wp,-2.03_wp,-2.05_wp,-2.06_wp,-2.08_wp, &
      -2.09_wp,-1.57_wp,-1.57_wp,-1.55_wp,-1.59_wp,-2.01_wp,-2.11_wp,-2.05_wp,-2.1_wp,-2.11_wp,-1.09_wp, &
      -1.14_wp,-1.08_wp,-1.13_wp, &
      -2.45_wp,-2.49_wp,-2.70_wp,-2.34_wp,-2.36_wp,-2.44_wp,-2.26_wp]
   real(wp), parameter :: BETASARR(32) = [ &
      -9.34_wp,-10.3_wp,-11.7_wp,-11.5_wp,-0.114_wp,-0.0407_wp,0.3_wp,0.668_wp,0.832_wp,0.995_wp, &
       1.29_wp,1.58_wp,1.1_wp,1.25_wp,1.33_wp,2.12_wp,0.894_wp,1.58_wp,1.19_wp,1.64_wp,2.1_wp,-0.37_wp, &
      -0.195_wp,-0.336_wp,-0.109_wp, &
       0.125_wp,0.345_wp,2.180_wp,-0.243_wp,-0.113_wp,0.254_wp,-3.460_wp]
   real(wp), parameter :: ATSARR(32) = [ &
      0.172_wp,0.174_wp,0.173_wp,0.171_wp,0.169_wp,0.166_wp,0.164_wp,0.189_wp,0.188_wp,0.185_wp, &
      0.184_wp,0.183_wp,0.198_wp,0.197_wp,0.195_wp,0.193_wp,0.198_wp,0.197_wp,0.197_wp,0.198_wp,0.198_wp, &
      0.218_wp,0.216_wp,0.216_wp,0.211_wp, &
      0.191_wp,0.184_wp,0.198_wp,0.184_wp,0.182_wp,0.188_wp,0.216_wp]
   real(wp), parameter :: CSARR(32) = [ &
      1.02e-12_wp,1.09e-12_wp,1.27e-12_wp,1.33e-12_wp,1.26e-13_wp,1.27e-13_wp,1.0e-13_wp, &
      5.2e-14_wp,4.81e-14_wp,4.7e-14_wp,4.26e-14_wp,3.94e-14_wp,4.24e-14_wp,4.0e-14_wp, &
      4.05e-14_wp,3.2e-14_wp,4.95e-14_wp,3.69e-14_wp,4.37e-14_wp,3.63e-14_wp,3.13e-14_wp, &
      1.17e-13_wp,1.05e-13_wp,1.17e-13_wp,1.04e-13_wp, &
      1.84e-14_wp,1.78e-14_wp,7.29e-15_wp,3.18e-14_wp,3.03e-14_wp,2.24e-14_wp,3.16e-14_wp]
   real(wp), parameter :: BC5ARR(32) = [ &
      0._wp,1._wp,2._wp,3._wp,4._wp,5._wp,6._wp,0._wp,1._wp,2._wp,3._wp,4._wp,0._wp,1._wp,2._wp,3._wp, &
      0._wp,1._wp,2._wp,3._wp,4._wp,0._wp,1._wp,2._wp,3._wp, &
      0._wp,1._wp,2._wp,0._wp,0.5_wp,1._wp,0._wp]

contains

   !====================================================================
   ! grain_dist_dl07: (1/n_H) dn/da [cm^-1] at radius a_um (microns).
   ! Faithful port of Qext_new/grain_dist.f90::GRAIN_DIST_DL07.
   !====================================================================
   function grain_dist_dl07(index, dtype, a_um) result(dnda)
      integer,          intent(in) :: index
      character(len=*), intent(in) :: dtype
      real(wp),         intent(in) :: a_um
      real(wp) :: dnda

      real(wp) :: a, alphag, betag, atg, acg, cg
      real(wp) :: alphas, betas, ats, acs, cs, bc5
      real(wp) :: dndavsg, a01, a02, sig1, sig2, b1, b2, xi_gra
      character(len=3) :: dt3
      integer :: i

      dnda = 0.0_wp
      if (index < 1 .or. index > 32) then
         write(*,*) 'grain_dist_dl07: index out of range [1,32]: ', index
         return
      end if

      ! Case-insensitive 3-char type
      dt3 = dtype(1:min(3,len(dtype)))
      do i = 1, 3
         if (dt3(i:i) >= 'A' .and. dt3(i:i) <= 'Z') &
            dt3(i:i) = achar(iachar(dt3(i:i)) + 32)
      end do

      a      = a_um * 1.0e-4_wp        ! micron -> cm
      alphag = ALPHAGARR(index)
      betag  = BETAGARR(index)
      atg    = ATGARR(index) * 1.0e-4_wp
      acg    = ACGARR(index) * 1.0e-4_wp
      cg     = CGARR(index)
      alphas = ALPHASARR(index)
      betas  = BETASARR(index)
      ats    = ATSARR(index) * 1.0e-4_wp
      acs    = 1.0e-5_wp
      cs     = CSARR(index)
      bc5    = BC5ARR(index)

      if (dt3 == 'sil') then
         dnda = (cs/a) * (a/ats)**alphas
         if (betas >= 0.0_wp) then
            dnda = dnda * (1.0_wp + betas*a/ats)
         else
            dnda = dnda / (1.0_wp - betas*a/ats)
         end if
         if (a > ats) dnda = dnda * exp(((ats-a)/acs)**3)

      else if (dt3 == 'gra' .or. dt3 == 'pah') then
         dnda = (cg/a) * (a/atg)**alphag
         if (betag >= 0.0_wp) then
            dnda = dnda * (1.0_wp + betag*a/atg)
         else
            dnda = dnda / (1.0_wp - betag*a/atg)
         end if
         if (a > atg) dnda = dnda * exp(((atg-a)/acg)**3)

         if (index <= 25) then
            ! Draine & Li (2007) VSG parameters (MW): matches DL07spec MW
            ! header (a_01=4e-8/sig=0.4, a_02=2e-7/sig=0.55).
            a01 = 4.0e-8_wp;  a02 = 2.0e-7_wp
            sig1 = 0.4_wp;    sig2 = 0.55_wp
            b1 = 1.2961e-7_wp; b2 = 1.2410e-10_wp
         else
            ! DL07 LMC2/SMC VSG parameters: DL07spec LMC2 (model 10) and SMC
            ! (model 11) headers both give a_01=4e-8/sig1=0.4, a_02=2.5e-7/
            ! sig2=0.4, with a 50/50 b_C split (b_C1=b_C2). The amplitudes
            ! b1,b2 are the WD01 eq.(4) lognormal normalizations for these
            ! (a0,sig) at b_C1=b_C2=5e-6, BC5=1 (validated vs the MW b1,b2 to
            ! 0.6%). [Previously used the older LD01 values 3.5e-8/3.0e-7,
            ! 2.0496e-7/9.6e-11, inconsistent with DL07.]
            a01 = 4.0e-8_wp;  a02 = 2.5e-7_wp
            sig1 = 0.4_wp;    sig2 = 0.4_wp
            b1 = 8.641e-8_wp; b2 = 3.318e-10_wp
         end if
         dndavsg = (b1/a) * exp(-0.5_wp*(log(a/a01)/sig1)**2) + &
                   (b2/a) * exp(-0.5_wp*(log(a/a02)/sig2)**2)
         if (dndavsg >= 1.0e-4_wp*dnda) dnda = dnda + bc5*dndavsg

         ! Split carbonaceous into graphite (xi_gra) and PAH (1-xi_gra)
         if (a <= 5.0e-7_wp) then
            xi_gra = 0.01_wp
         else
            xi_gra = 0.01_wp + 0.99_wp*(1.0_wp - (5.0e-7_wp/a)**3)
         end if
         if (dt3 == 'gra') then
            dnda = xi_gra * dnda
         else
            dnda = (1.0_wp - xi_gra) * dnda
         end if

      else
         write(*,*) 'grain_dist_dl07: dtype should be sil, gra, or pah; got ', trim(dtype)
      end if

      ! Draine (2003, ARA&A) abundance reduction for the MW sets. Disabled
      ! by default (gd_apply_d03_reduction=.false.) to match the
      ! DL07 published spectra, which use the un-reduced WD01 abundances.
      if (index <= 25 .and. gd_apply_d03_reduction) dnda = dnda * 0.93_wp
   end function grain_dist_dl07

end module grain_dist_mod
