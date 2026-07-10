module mc_grain_type
   ! Thread-safe container for single-grain state of the MC P(T) simulation.
   ! Holds the pre-tabulated absorption cross section, the local radiation
   ! field, the photon-energy sampling CDF, and the temperature-dependent
   ! tables (enthalpy, Planck-averaged Q).  One mc_grain_t per grain;
   ! distinct instances may live on distinct OpenMP threads without
   ! locking.

   use constants, only: wp
   implicit none
   private
   public :: mc_grain_t

   type :: mc_grain_t
      ! ---- grain identity ----
      real(wp)          :: a_um   = -1.0_wp       ! grain radius [um]
      character(len=16) :: comp   = ''            ! heat-capacity backend tag

      ! ---- wavelength grid + Cabs / u_lam ----
      integer  :: NLAM = 0
      real(wp), allocatable :: lam_grid(:)        ! [um], log-spaced
      real(wp), allocatable :: Q_grid(:)          ! Q_abs(a, lam)
      real(wp), allocatable :: u_lam_grid(:)      ! [erg / cm^3 / um]

      ! ---- continuous + stochastic split ----
      real(wp) :: lam_c     = 1000.0_wp           ! [um]
      real(wp) :: H_cont    = 0.0_wp              ! [erg / cm^2 / s]
      real(wp) :: rate_event= 0.0_wp              ! [1/s]

      ! ---- photon-wavelength sampling CDF ----
      integer  :: NCDF = 0
      real(wp), allocatable :: cdf_F(:)
      real(wp), allocatable :: cdf_lam(:)

      ! ---- thermo tables ----
      integer  :: NT = 0
      real(wp), allocatable :: T_grid(:)          ! [K], log-spaced
      real(wp), allocatable :: U_grid(:)          ! enthalpy [erg]
      real(wp), allocatable :: QPL_grid(:)        ! <Q>_T
   end type mc_grain_t

end module mc_grain_type
