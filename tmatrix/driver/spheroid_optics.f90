module spheroid_optics
   ! Random-orientation optical efficiencies (Q_ext, Q_abs, Q_sca, albedo, g)
   ! of a homogeneous spheroid over the whole size-parameter range.
   !
   ! The physics owned here is the choice of regime.  With
   !
   !     x = 2 pi a_eff / lambda
   !
   ! the response is evaluated in one of three domains of validity:
   !
   !   x < X_RAYLEIGH_MAX (0.1)      Rayleigh dipole limit.  The grain is much
   !                                 smaller than the wavelength, the internal
   !                                 field is uniform, and the orientation-
   !                                 averaged polarizability of the spheroid
   !                                 gives Q_abs ~ x and Q_sca ~ x^4 in closed
   !                                 form.  The T-matrix expansion is also
   !                                 poorly conditioned here: the scattered
   !                                 field falls far below the convergence
   !                                 tolerance of the field expansion.
   !
   !   X_RAYLEIGH_MAX <= x <= X_GEOMETRIC_MIN (50)
   !                                 T-matrix.  Waterman's extended boundary
   !                                 condition solution for a spheroid with the
   !                                 orientation average taken analytically.
   !                                 This is the only regime that resolves the
   !                                 interference and ripple structure of Q(x),
   !                                 which neither asymptotic form reproduces.
   !
   !   x > X_GEOMETRIC_MIN           Geometric optics limit.  The grain is much
   !                                 larger than the wavelength: Q_ext -> 2
   !                                 (extinction paradox) and absorption follows
   !                                 the attenuation along a single chord.  The
   !                                 required T-matrix expansion order grows
   !                                 roughly like x, so both the storage bounds
   !                                 and the conditioning of the solver become
   !                                 limiting around here.
   !
   ! When the T-matrix solver reports a nonzero IERR the result is taken from
   ! whichever asymptotic limit the size parameter is closer to, split at
   ! X_LIMIT_SPLIT (x = 1, where the grain circumference equals the wavelength).
   !
   ! flag convention returned by spheroid_q (unchanged from the original
   ! astrodust table driver, so previously written tables stay readable):
   !
   !    0     T-matrix converged
   !   10     Rayleigh dipole limit
   !   20     geometric optics limit
   !   11..19 T-matrix returned IERR = 1..9 and x <  X_LIMIT_SPLIT, so the row
   !          was taken from the Rayleigh dipole limit (flag = IERR + 10)
   !   21..29 T-matrix returned IERR = 1..9 and x >= X_LIMIT_SPLIT, so the row
   !          was taken from the geometric optics limit (flag = IERR + 20)
   !
   ! A caller that also runs check_physical_bounds is expected to add 100 to the
   ! flag of a row that fails it, preserving the base flag in the low two digits.
   !
   ! Conventions.  Refractive index is m = nr + i*ki with ki >= 0; a_um and
   ! lambda_um share the same length unit (micron in the astrodust tables), and
   ! the returned efficiencies are normalized to pi*a_eff^2.  Shape follows
   ! Mishchenko: options%shape = -1 is a spheroid and
   ! options%aspect_ratio = b/a (equatorial axis over rotational axis), so
   ! aspect_ratio > 1 is oblate, < 1 prolate, = 1 a sphere.
   !
   ! This module depends only on the T-matrix library and on the closed-form
   ! asymptotic limits.  It reads no dielectric table and knows nothing about
   ! any particular grain material or grid.
   !
   ! geometric_optics_limit is called here with the axis ratio eps_ba =
   ! options%aspect_ratio.  The averaged Q_ext = 2, Q_abs = 1 - exp(-4 k x),
   ! Q_sca and g = 0 returned to spheroid_q do not contain eps_ba at all; the
   ! argument is consumed only by that routine's optional orientation-resolved
   ! outputs, which spheroid_q does not request.  Passing it therefore changes
   ! no number here.

   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use tmatrix_kinds, only: wp
   use tmatrix_status, only: TMATRIX_SUCCESS
   use tmatrix_api, only: tmatrix_options_t, tmatrix_result_t, &
                          tmatrix_workspace_t, tmatrix_eval
   use asymptotic_optics, only: rayleigh_limit, geometric_optics_limit
   implicit none
   private

   public :: spheroid_q, check_physical_bounds
   public :: X_RAYLEIGH_MAX, X_GEOMETRIC_MIN, X_LIMIT_SPLIT

   !! Upper bound of the Rayleigh dipole approximation.  The electric-dipole
   !! term dominates and the neglected multipoles enter at O(x^2), the
   !! conventional few-percent boundary for |m| x << 1.
   real(wp), parameter :: X_RAYLEIGH_MAX  = 0.1_wp
   !! Lower bound of the geometric optics approximation, where the grain is
   !! opaque on the scale of the wavelength and Q_ext has reached its
   !! diffraction-plus-interception value of 2.
   real(wp), parameter :: X_GEOMETRIC_MIN = 50.0_wp
   !! Size parameter at which the grain circumference equals the wavelength.
   !! It separates the two asymptotic limits when the T-matrix solver fails.
   real(wp), parameter :: X_LIMIT_SPLIT   = 1.0_wp

   real(wp), parameter :: PI = acos(-1.0_wp)

contains

   subroutine spheroid_q(work, a_um, lambda_um, nr, ki, options, &
                         qext, qabs, qsca, albedo, asymmetry, flag, status, message)
      !! Orientation-averaged efficiencies of one spheroid at one wavelength.
      !! The regime is selected from the size parameter as described in the
      !! module header; flag records which one produced the answer.
      !!
      !! status and message pass the T-matrix library verdict through unchanged.
      !! A nonzero status with flag = 0 means the library reported a failure the
      !! IERR redirection does not cover (invalid input, uninitialized or
      !! unallocatable workspace, or a nonfinite result); the returned values are
      !! then whatever the library produced and the caller should submit them to
      !! check_physical_bounds before use.
      type(tmatrix_workspace_t), intent(inout) :: work
      real(wp), intent(in)  :: a_um, lambda_um, nr, ki
      type(tmatrix_options_t), intent(in) :: options
      real(wp), intent(out) :: qext, qabs, qsca, albedo, asymmetry
      integer,  intent(out) :: flag
      integer,          intent(out), optional :: status
      character(len=*), intent(out), optional :: message

      type(tmatrix_result_t) :: result
      real(wp) :: x
      integer  :: ierr_t

      x = 2.0_wp * PI * a_um / lambda_um

      if (x < X_RAYLEIGH_MAX) then
         call rayleigh_limit(a_um, lambda_um, nr, ki, options%aspect_ratio, &
                             qext, qsca, albedo, asymmetry)
         qabs = qext - qsca
         flag = 10
         if (present(status))  status  = TMATRIX_SUCCESS
         if (present(message)) message = ''
      else if (x > X_GEOMETRIC_MIN) then
         call geometric_optics_limit(a_um, lambda_um, nr, ki, options%aspect_ratio, &
                                     qext, qsca, albedo, asymmetry)
         qabs = qext - qsca
         flag = 20
         if (present(status))  status  = TMATRIX_SUCCESS
         if (present(message)) message = ''
      else
         call tmatrix_eval(work, a_um, lambda_um, nr, ki, options, result)
         qext      = result%qext
         qsca      = result%qsca
         albedo    = result%albedo
         asymmetry = result%asymmetry
         ierr_t    = result%legacy_ierr
         if (ierr_t /= 0) then
            if (x < X_LIMIT_SPLIT) then
               call rayleigh_limit(a_um, lambda_um, nr, ki, options%aspect_ratio, &
                                   qext, qsca, albedo, asymmetry)
               flag = ierr_t + 10
            else
               call geometric_optics_limit(a_um, lambda_um, nr, ki, options%aspect_ratio, &
                                           qext, qsca, albedo, asymmetry)
               flag = ierr_t + 20
            end if
         else
            flag = 0
         end if
         qabs = qext - qsca
         if (present(status))  status  = result%status
         if (present(message)) message = trim(result%message)
      end if
   end subroutine spheroid_q


   subroutine check_physical_bounds(qext, qabs, qsca, albedo, asymmetry, ok, reason)
      ! Physical-consistency test for one (lambda, a_eff) result.
      ! Sets ok = .false. and appends a short tag to reason for every bound
      ! that is violated; reason is empty when everything passes.
      real(wp),         intent(in)  :: qext, qabs, qsca, albedo, asymmetry
      logical,          intent(out) :: ok
      character(len=*), intent(out) :: reason
      real(wp), parameter :: tol = 1.0e-9_wp

      reason = ''
      ! Finiteness (NaN / Inf).
      if (.not. ieee_is_finite(qext))      reason = trim(reason)//' qext-nonfinite'
      if (.not. ieee_is_finite(qsca))      reason = trim(reason)//' qsca-nonfinite'
      if (.not. ieee_is_finite(qabs))      reason = trim(reason)//' qabs-nonfinite'
      if (.not. ieee_is_finite(albedo))    reason = trim(reason)//' albedo-nonfinite'
      if (.not. ieee_is_finite(asymmetry)) reason = trim(reason)//' g-nonfinite'
      ! Non-negativity (small negative roundoff tolerated) and qext >= qsca.
      if (ieee_is_finite(qext) .and. qext < -tol) reason = trim(reason)//' qext<0'
      if (ieee_is_finite(qsca) .and. qsca < -tol) reason = trim(reason)//' qsca<0'
      if (ieee_is_finite(qabs) .and. qabs < -tol) reason = trim(reason)//' qabs<0'
      if (ieee_is_finite(qext) .and. ieee_is_finite(qsca) .and. qext < qsca - tol) &
         reason = trim(reason)//' qext<qsca'
      ! Albedo in [0,1] and asymmetry g in [-1,1].
      if (ieee_is_finite(albedo) .and. (albedo < -tol .or. albedo > 1.0_wp + tol)) &
         reason = trim(reason)//' albedo-range'
      if (ieee_is_finite(asymmetry) .and. &
          (asymmetry < -1.0_wp - tol .or. asymmetry > 1.0_wp + tol)) &
         reason = trim(reason)//' g-range'

      ok = (len_trim(reason) == 0)
   end subroutine check_physical_bounds

end module spheroid_optics
