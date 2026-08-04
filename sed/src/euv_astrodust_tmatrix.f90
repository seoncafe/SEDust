module euv_astrodust_tmatrix
   ! Astrodust optics of the extreme-ultraviolet band, i.e. of the wavelengths
   ! sed_init prepends below the T-matrix Q table's 0.0912 um (13.6 eV) end.
   !
   ! The band is computed on the SAME particle the table is: the b/a = 1.400
   ! oblate spheroid, orientation-averaged, over the whole size-parameter
   ! range.  spheroid_q picks the regime from x = 2 pi a_eff / lambda -- the
   ! Rayleigh dipole limit below x = 0.1, the T-matrix between, the geometric
   ! optics limit above x = 50 -- exactly as the table's own driver does, so
   ! the band continues the table's calculation instead of gluing a different
   ! grain to its short-wavelength end.
   !
   ! This is the ONE file of the SED library that reaches into libtmatrix.a.
   ! sed_astrodust_mod names the calculation through a procedure pointer, not a
   ! use statement, so a library built without this file carries no T-matrix
   ! reference and links without libtmatrix.a; such a build refuses
   ! euv_tmatrix = .true. (status 6) rather than substituting the
   ! volume-equivalent sphere for the spheroid.  Compile this file, link
   ! libtmatrix.a, and call
   !
   !    call use_tmatrix_euv_band_optics()
   !
   ! once before build_astrodust to get the spheroid route.
   !
   ! The refractive index arrives as an argument, so nothing here reads a
   ! dielectric table or touches the state of the SED solver.
   use tmatrix_kinds,     only: wp
   use constants,         only: sed_wp => wp
   use tmatrix_api,       only: tmatrix_options_t, tmatrix_workspace_t, &
                                tmatrix_workspace_init, tmatrix_workspace_finalize
   use spheroid_optics,   only: spheroid_q, check_physical_bounds
   use sed_astrodust_mod, only: sed_register_euv_band_optics
   implicit none
   private
   public :: spheroid_euv_band_optics, use_tmatrix_euv_band_optics

   ! libtmatrix.a owns its working precision (tmatrix_kinds) and the SED side
   ! owns its own (constants), on purpose: neither imports the other's.  The
   ! routine below is declared with the T-matrix kind and called through an
   ! interface declared with the SED kind, so the two must be the same kind.
   ! They are; this division stops the build on the day they are not.
   integer, parameter :: WP_KINDS_AGREE = 1 / merge(1, 0, wp == sed_wp)

   ! Particle the astrodust EUV band is solved for.  These are the settings
   ! tmatrix/driver/run_tmatrix.f90 used to compute
   ! q_astrodust_P0.20_Fe0.00_1.400.dat (its EPS_BA, DDELT, NP_OBL, NDGS), and
   ! they must stay equal to them: the extension and the table are then one
   ! calculation of one particle, differing only in wavelength.  A mismatch
   ! would change the grain's shape or the solver's accuracy at the seam.
   real(wp), parameter :: AD_AXIAL_RATIO = 1.4_wp      ! b/a > 1: oblate spheroid
   real(wp), parameter :: AD_TM_TOL      = 1.0e-3_wp   ! T-matrix convergence tolerance
   integer,  parameter :: AD_TM_SHAPE    = -1          ! Mishchenko code for a spheroid
   integer,  parameter :: AD_TM_NDGS     = 2           ! Gaussian division points per order

contains

   subroutine use_tmatrix_euv_band_optics()
      !! Make the spheroid the grain sed_init solves the EUV band for.  One
      !! call, any time before build_astrodust / sed_init.
      call sed_register_euv_band_optics(spheroid_euv_band_optics)
   end subroutine use_tmatrix_euv_band_optics


   subroutine spheroid_euv_band_optics(lam, a_um, nr, ki, qabs, qsca, gsca, status)
      !! Orientation-averaged efficiencies of the oblate spheroid at every
      !! (wavelength, radius) pair of the EUV band.  Implements
      !! sed_astrodust_mod's euv_band_optics_i; see there for the arguments
      !! and the meaning of status.
      real(wp), intent(in)  :: lam(:)          ! (n_euv) [um]
      real(wp), intent(in)  :: a_um(:)         ! (n_size) [um]
      real(wp), intent(in)  :: nr(:), ki(:)    ! (n_euv) m = nr + i*ki
      real(wp), intent(out) :: qabs(:,:)       ! (n_euv, n_size)
      real(wp), intent(out) :: qsca(:,:)
      real(wp), intent(out) :: gsca(:,:)
      integer,  intent(out) :: status
      integer  :: n_euv, n_size, ja, jw, n_bad
      integer  :: flag, tm_stat
      real(wp) :: qext1, qsca1, qabs1, alb1, gsca1
      logical  :: phys_ok
      character(len=160) :: tm_msg
      character(len=192) :: phys_why
      type(tmatrix_options_t)   :: tm_opts
      type(tmatrix_workspace_t) :: tm_work

      status = 0
      n_euv  = size(lam)
      n_size = size(a_um)
      if (size(nr) /= n_euv .or. size(ki) /= n_euv .or.               &
          size(qabs, 1) /= n_euv .or. size(qabs, 2) /= n_size .or.    &
          size(qsca, 1) /= n_euv .or. size(qsca, 2) /= n_size .or.    &
          size(gsca, 1) /= n_euv .or. size(gsca, 2) /= n_size) then
         status = -1
         return
      end if
      if (n_euv <= 0 .or. n_size <= 0) return

      tm_opts%aspect_ratio = AD_AXIAL_RATIO
      tm_opts%tolerance    = AD_TM_TOL
      tm_opts%shape        = AD_TM_SHAPE
      tm_opts%ndgs         = AD_TM_NDGS
      n_bad = 0

      ! Threading.  Each ja writes its own column of qabs / qsca / gsca and
      ! reads nothing another ja writes, so the loop carries no dependence and
      ! the result cannot depend on the schedule or the thread count.
      !$omp parallel default(none) &
      !$omp&   shared(n_euv, n_size, lam, a_um, nr, ki, qabs, qsca, gsca, tm_opts) &
      !$omp&   private(ja, jw, qext1, qsca1, qabs1, alb1, gsca1, &
      !$omp&           flag, tm_stat, tm_msg, tm_work, phys_ok, phys_why) &
      !$omp&   reduction(+:n_bad)
      ! One T-matrix workspace for each thread, because a workspace may have
      ! only one active caller.  It costs 37 MiB here and grows to 80 MiB on
      ! the first evaluation, so nothing is allocated until a band is asked
      ! for.
      call tmatrix_workspace_init(tm_work)
      !$omp do schedule(dynamic)
      do ja = 1, n_size
         do jw = 1, n_euv
            call spheroid_q(tm_work, a_um(ja), lam(jw), nr(jw), ki(jw), tm_opts, &
                            qext1, qabs1, qsca1, alb1, gsca1, flag, &
                            status=tm_stat, message=tm_msg)
            call check_physical_bounds(qext1, qabs1, qsca1, alb1, gsca1, &
                                       phys_ok, phys_why)
            ! A nonzero status with flag = 0 is a library failure the
            ! asymptotic-limit redirection does not cover, so report it too.
            if (tm_stat /= 0 .and. flag == 0) then
               phys_ok  = .false.
               phys_why = trim(phys_why)//' '//trim(tm_msg)
            end if
            if (.not. phys_ok) then
               ! Warn and carry on: this may be running inside an RT host,
               ! which is told about failures through `status`, not by having
               ! its process stopped.
               n_bad = n_bad + 1
               !$omp critical (euv_optics_warning)
               write(*,'(a,es12.5,a,es12.5,a,i0,a,a)') &
                  ' euv_astrodust_tmatrix: EUV optics failed a physical bound at lambda=', &
                  lam(jw), ' um, a_eff=', a_um(ja), ' um, flag=', flag, &
                  ', broken:', trim(phys_why)
               !$omp end critical (euv_optics_warning)
            end if
            qabs(jw, ja) = qabs1
            qsca(jw, ja) = qsca1
            gsca(jw, ja) = gsca1
         end do
      end do
      !$omp end do
      call tmatrix_workspace_finalize(tm_work)
      !$omp end parallel

      status = n_bad
   end subroutine spheroid_euv_band_optics

end module euv_astrodust_tmatrix
