module constants
public

! kind parameter for single and double precision (6 and 15 significant decimal digits)
! precision of 32-, 64-, and 128-bit reals
  integer, parameter :: sp = selected_real_kind(6, 37)
  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: qp = selected_real_kind(33, 4931)

! precision of 32-, 64- bit integers
  integer, parameter :: i4b = selected_int_kind(r=9)
  integer, parameter :: i8b = selected_int_kind(r=18)

! working precision
! Note that (1) iwp (integer working precision) applies only to random number seed.
!           (2) nphotons is always 64-bit integer.
!           (3) other integer variables are always 32-bit integers.
!  integer, parameter :: wp  = sp
!  integer, parameter :: iwp = i4b
  integer, parameter :: wp  = dp
  integer, parameter :: iwp = i8b

! tinest = the smallest positive number
! eps    = the least positive number
!          that added to 1 returns a number that is greater than 1
! hugest = the hugest positive number
  real(kind=wp), parameter :: tinest = tiny(0.0_wp)
  real(kind=wp), parameter :: eps    = epsilon(0.0_wp)
  real(kind=wp), parameter :: hugest = huge(1.0_wp)

! numerical values
  real(kind=wp), parameter :: pi     = 3.141592653589793238462643383279502884197_wp
  real(kind=wp), parameter :: twopi  = 6.283185307179586476925286766559005768394_wp
  real(kind=wp), parameter :: fourpi = 12.56637061435917295385057353311801153679_wp

  real(kind=wp), parameter :: onethird = 1.0_wp/3.0_wp
  real(kind=wp), parameter :: twothird = 2.0_wp/3.0_wp

! conversion factor (radian to degree, degree to radian)
  real(kind=wp), parameter :: rad2deg = 180.0_wp/pi
  real(kind=wp), parameter :: deg2rad = pi/180.0_wp

! distances
  real(kind=wp), parameter :: um2cm  = 1e-4_wp
  real(kind=wp), parameter :: cm2um  = 1e4_wp
  real(kind=wp), parameter :: pc2cm  = 3.0856776e18_wp
  real(kind=wp), parameter :: kpc2cm = pc2cm * 1e3_wp
  real(kind=wp), parameter :: au2cm  = 1.4960e13_wp

! physical constants, cgs.
! h, c, k_B and e are the defining constants of the SI (exact since the 2019
! redefinition; BIPM, The International System of Units, 9th ed., 2019, Table 1).
! Everything below them is derived from them here, so no module carries its own
! rounded copy. The SI-unit values are the same exact numbers.
  real(kind=wp), parameter :: c        = 2.99792458e10_wp     ! speed of light [cm s^-1]
  real(kind=wp), parameter :: h_planck = 6.62607015e-27_wp    ! Planck constant [erg s]
  real(kind=wp), parameter :: k_boltz  = 1.380649e-16_wp      ! Boltzmann constant [erg K^-1]
  real(kind=wp), parameter :: ev2erg   = 1.602176634e-12_wp   ! 1 eV [erg] (e = 1.602176634e-19 C)

  real(kind=wp), parameter :: c_si     = 2.99792458e8_wp      ! [m s^-1]
  real(kind=wp), parameter :: h_si     = 6.62607015e-34_wp    ! [J s]
  real(kind=wp), parameter :: k_si     = 1.380649e-23_wp      ! [J K^-1]

  real(kind=wp), parameter :: hc       = h_planck * c         ! [erg cm]; also erg per cm^-1 of wavenumber
  real(kind=wp), parameter :: hc2      = hc**2                ! [erg^2 cm^2]
  real(kind=wp), parameter :: hc_ergum = hc * 1e4_wp          ! [erg um]
  real(kind=wp), parameter :: hc_evum  = hc_ergum / ev2erg    ! [eV um] = 1.2398419843...
  real(kind=wp), parameter :: hc_k     = hc / k_boltz         ! second radiation constant [cm K] = 1.4387768775...

! dust mass per hydrogen nuclei (g/H)
  real(kind=wp), parameter :: h2dust = 1.87e-26_wp

! kappa_v : dust extinction cross-section (cm^2/g) at V-band
  real(kind=wp), parameter :: kappa_v = 25992.6_wp

end module constants
