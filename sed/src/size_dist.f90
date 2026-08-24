module size_dist_mod
   ! HD23 (Hensley & Draine 2023, ApJ 948, 55) astrodust+PAH size
   ! distributions, evaluated from their analytic forms:
   !   eq. (24)  dn_Ad/da   lognormal + quintic log-polynomial, 4.5 A < a < 5 um
   !   eq. (17)  dn_PAH/da  two lognormals, a > 4.0 A
   !   eq. (19)  f_ion(a)   PAH ionized fraction
   !   eq. (25)  f_align(a) astrodust alignment fraction
   ! on the 167-point grid of the HD23 release table
   ! (log10 a/um = -3.45 .. 0.70, step 0.025), with dn_i/nH already
   ! integrated over each bin, da = a (10^0.0125 - 10^-0.0125).  Do not
   ! multiply by da again.
   !
   ! CONSTANTS.  The values printed in HD23 Table 1 carry three significant
   ! figures.  That is not enough for eq. (24): the quintic term
   ! A_5 (ln a)^5 reaches -247 at a = 5 um, so the rounding of A_5 alone
   ! changes dn_Ad/da there by a factor ~2.  Evaluated with the printed
   ! values, dn_Ad exceeds the release table by 2% at 0.1 um, 9% at 0.3 um,
   ! 81% at 5 um, and by 20% in total volume.  The constants below are the
   ! same functional forms refit to the release table
   ! (data/release/size_distribution.dat).  They agree with Table 1 at its
   ! three figures and reproduce the table to |dlog| <= 1.6e-4, the table's
   ! own five-digit rounding.  The PAH constants a_0,j, sigma and the
   ! ionization law are exact as printed; B_2 and the alignment constants
   ! were refit the same way (their printed rounding costs 2e-4 and 3e-3).
   !
   ! test_size_distribution.x compares this routine against the table.

   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: hd23_size_distribution
   public :: n_size, a_dist, dn_ad, dn_pah, f_ion, f_align
   ! The HD23 alignment law, eq. (25), for any radius; and the same power
   ! law with its three parameters free, for a host imposing its own state.
   public :: falign_hd23, falign_powerlaw
   public :: A_ALIGN, ALPHA_ALIGN, FMAX_ALIGN

   integer, parameter :: wp = real64

   integer  :: n_size = 0
   real(wp), allocatable :: a_dist(:)   ! [microns], ascending
   real(wp), allocatable :: dn_ad(:)    ! [1/H]  (already binned)
   real(wp), allocatable :: dn_pah(:)   ! [1/H]
   real(wp), allocatable :: f_ion(:)    ! PAH ionization fraction
   real(wp), allocatable :: f_align(:)  ! alignment fraction

   ! Grid of the HD23 release table.
   real(wp), parameter :: LOG10_AMIN = -3.45_wp, DLOG10 = 0.025_wp
   integer,  parameter :: NPT = 167

   ! eq. (24), astrodust.  Radii in Angstrom.
   real(wp), parameter :: A_AD_MIN_ANG = 4.5_wp, A_AD_MAX_UM = 5.0_wp
   real(wp), parameter :: B_AD   = 3.311789e-10_wp        ! Table 1: 3.31e-10 H^-1
   real(wp), parameter :: A0_AD  = 63.8091_wp             !          63.8 A
   real(wp), parameter :: SIG_AD = 0.352558_wp            !          0.353
   real(wp), parameter :: A_POLY(0:5) = [ 2.9747134e-5_wp, &   !  2.97e-5 H^-1
                                         -3.4022298_wp,   &   ! -3.40
                                         -0.80690301_wp,  &   ! -0.807
                                          0.15654565_wp,  &   !  0.157
                                          7.9647975e-3_wp, &  !  7.96e-3
                                         -1.6804908e-3_wp ]   ! -1.68e-3
   ! eq. (17), PAH.
   real(wp), parameter :: A_PAH_MIN_ANG = 4.0_wp
   real(wp), parameter :: A0_PAH(2)  = [4.0_wp, 30.0_wp]
   real(wp), parameter :: SIG_PAH    = 0.40_wp
   real(wp), parameter :: B_PAH(2)   = [7.52e-7_wp, 8.09137799e-10_wp]  ! 7.52e-7, 8.09e-10
   ! eq. (19), PAH ionized fraction.
   real(wp), parameter :: A_ION_ANG = 10.0_wp
   ! eq. (25), alignment:  f_align(a) = f_max / (1 + (a_align/a)**alpha_align)
   real(wp), parameter :: A_ALIGN     = 0.074920_wp        ! Table 1: 0.0749 um
   real(wp), parameter :: ALPHA_ALIGN = 1.799319_wp        !          1.80
   real(wp), parameter :: FMAX_ALIGN  = 1.0_wp             !          1.00

contains

   subroutine hd23_size_distribution()
      real(wp) :: a_um, a_ang, da_ang, lna, lognormal, poly
      integer  :: i, k

      n_size = NPT
      if (allocated(a_dist))  deallocate(a_dist)
      if (allocated(dn_ad))   deallocate(dn_ad)
      if (allocated(dn_pah))  deallocate(dn_pah)
      if (allocated(f_ion))   deallocate(f_ion)
      if (allocated(f_align)) deallocate(f_align)
      allocate(a_dist(NPT), dn_ad(NPT), dn_pah(NPT), f_ion(NPT), f_align(NPT))

      do i = 1, NPT
         a_um   = 10.0_wp**(LOG10_AMIN + DLOG10*(i-1))
         a_ang  = a_um * 1.0e4_wp
         da_ang = a_ang * (10.0_wp**(DLOG10/2) - 10.0_wp**(-DLOG10/2))
         lna    = log(a_ang)
         a_dist(i) = a_um

         if (a_ang > A_AD_MIN_ANG .and. a_um < A_AD_MAX_UM) then
            lognormal = B_AD/a_ang * exp(-(log(a_ang/A0_AD))**2 / (2*SIG_AD**2))
            poly = 0.0_wp
            do k = 1, 5
               poly = poly + A_POLY(k)*lna**k
            end do
            dn_ad(i) = (lognormal + A_POLY(0)/a_ang * exp(poly)) * da_ang
         else
            dn_ad(i) = 0.0_wp
         end if

         if (a_ang > A_PAH_MIN_ANG) then
            dn_pah(i) = 0.0_wp
            do k = 1, 2
               dn_pah(i) = dn_pah(i) + B_PAH(k)/a_ang * &
                  exp(-(log(a_ang/A0_PAH(k)))**2 / (2*SIG_PAH**2))
            end do
            dn_pah(i) = dn_pah(i) * da_ang
         else
            dn_pah(i) = 0.0_wp
         end if

         f_ion(i)   = 1.0_wp - 1.0_wp/(1.0_wp + a_ang/A_ION_ANG)
         f_align(i) = falign_hd23(a_um)
      end do
   end subroutine hd23_size_distribution


   ! HD23 eq. (25) for a grain of effective radius a [um].
   pure function falign_hd23(a) result(f)
      real(wp), intent(in) :: a
      real(wp)             :: f
      f = falign_powerlaw(a, FMAX_ALIGN, A_ALIGN, ALPHA_ALIGN)
   end function falign_hd23


   ! The same power-law rolloff with the three parameters left free,
   !
   !   f_align(a) = f_max / (1 + (a_align/a)**alpha_align),
   !
   ! so that a host can impose an alignment state other than the HD23 fit --
   ! a different critical radius a_align where the efficiency reaches
   ! f_max/2, a different sharpness alpha_align, or an overall reduction
   ! f_max < 1 standing for imperfect alignment. a is the effective radius
   ! [um]; f = 0 for a non-positive radius.
   pure function falign_powerlaw(a, f_max, a_align, alpha_align) result(f)
      real(wp), intent(in) :: a, f_max, a_align, alpha_align
      real(wp)             :: f
      if (a <= 0.0_wp) then
         f = 0.0_wp
      else
         f = f_max / (1.0_wp + (a_align / a)**alpha_align)
      end if
   end function falign_powerlaw

end module size_dist_mod
