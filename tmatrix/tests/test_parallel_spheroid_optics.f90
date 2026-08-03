program test_parallel_spheroid_optics
   !! Driver-layer concurrency test.  The other parallel harnesses stop at the
   !! library boundary (`tmatrix_eval`, `ampl`); this one enters `spheroid_q`,
   !! which is what a host actually calls from an OpenMP loop, and through it
   !! the two closed-form limits of `asymptotic_optics`.
   !!
   !! Coverage is by size-parameter regime, x = 2 pi a_eff / lambda:
   !!
   !!   x <  0.1   Rayleigh dipole limit      -> asymptotic_optics
   !!   0.1 <= x <= 50   T-matrix             -> tmatrix_eval
   !!   x >  50    geometric optics limit     -> asymptotic_optics
   !!
   !! The two asymptotic branches are the ones that leave the library entirely,
   !! so they are the reason this harness exists.  Each is additionally called
   !! on its own, with the orientation-resolved outputs requested, because
   !! those arguments are what reach `geometric_optics_limit`'s local arrays
   !! and the surface quadrature of `fresnel_opaque_absorption` -- storage that
   !! `spheroid_q` itself never activates.  Those direct calls evaluate the
   !! closed forms at every case, including sizes outside the domain of
   !! validity of the limit being evaluated: the quantity under test is
   !! reproducibility of the arithmetic, not the physics of the result.
   !!
   !! Every routine here is deterministic and reads only its own arguments and
   !! its caller-owned workspace, so serial and concurrent runs must agree in
   !! every bit.  The comparison is on IEEE bit patterns, not a tolerance: a
   !! difference means one thread's locals were another thread's locals.
   use, intrinsic :: iso_fortran_env, only: real64, int64
   use omp_lib, only: omp_get_thread_num, omp_get_max_threads
   use tmatrix_api, only: tmatrix_options_t, tmatrix_workspace_t, &
                          tmatrix_workspace_init, tmatrix_workspace_finalize
   use spheroid_optics, only: spheroid_q
   use asymptotic_optics, only: rayleigh_limit, geometric_optics_limit
   implicit none

   !! Geometry catalogue.  Four rows fall in the Rayleigh regime, six in the
   !! T-matrix regime, and two in the geometric-optics regime; the T-matrix
   !! rows stay inside the storage limits of the solver so the test measures
   !! reentrancy rather than overflow behavior.
   integer, parameter :: n_geom = 12
   real(real64), parameter :: a_list(n_geom) = &
      [0.0004_real64, 0.0020_real64, 0.0040_real64, 0.0080_real64, &
       0.0080_real64, 0.0500_real64, 0.0500_real64, 0.2500_real64, &
       0.2500_real64, 0.5000_real64, 2.5000_real64, 10.000_real64]
   real(real64), parameter :: lam_list(n_geom) = &
      [1.0000_real64, 1.0000_real64, 0.5000_real64, 0.5500_real64, &
       0.2500_real64, 0.5500_real64, 0.2500_real64, 0.5500_real64, &
       0.2500_real64, 0.2500_real64, 0.2500_real64, 0.2500_real64]
   !! x = 2 pi a / lambda for the rows above:
   !!   0.0025 0.0126 0.0503 0.0914 | 0.201 0.571 1.257 2.856 6.283 12.566 |
   !!   62.83 251.3
   real(real64), parameter :: absorption(2) = [0.001_real64, 0.100_real64]
   real(real64), parameter :: aspect(2)     = [1.400_real64, 0.700_real64]
   integer,      parameter :: n_cases = n_geom * 4
   integer,      parameter :: spheroid_shape = -1
   real(real64), parameter :: tolerance = 1.0e-3_real64
   integer,      parameter :: ndgs = 2
   !! Length of the coefficient arrays the analytic expansions write (LMAX = 2
   !! for the dipole matrix, 0 for the isotropic geometric-optics matrix).
   integer, parameter :: n_coeff = 3
   integer, parameter :: n_ori = 3

   type :: case_record_t
      !! Everything one case produces, in one place, so serial and parallel
      !! runs are compared field by field.
      real(real64) :: qext, qabs, qsca, albedo, asymmetry
      integer      :: flag, status
      real(real64) :: ray_q(4), go_q(4)
      real(real64) :: ray_coeff(n_coeff,6), go_coeff(n_coeff,6)
      integer      :: ray_lmax, go_lmax
      real(real64) :: ray_ori(n_ori,4), go_ori(n_ori,4)
   end type case_record_t

   type(tmatrix_workspace_t), allocatable :: work(:)
   type(case_record_t) :: serial_record(n_cases), parallel_record(n_cases)
   integer :: i, nwork, failures
   integer :: n_rayleigh, n_tmatrix, n_geometric

   nwork = max(1, omp_get_max_threads())
   allocate(work(nwork))
   do i = 1, nwork
      call tmatrix_workspace_init(work(i))
      if (.not. work(i)%initialized) error stop 'spheroid-optics workspace allocation failed'
   end do

   do i = 1, n_cases
      call evaluate(work(1), i, serial_record(i))
   end do

   !$omp parallel do default(none) private(i) shared(work,parallel_record) schedule(static)
   do i = 1, n_cases
      call evaluate(work(omp_get_thread_num()+1), i, parallel_record(i))
   end do
   !$omp end parallel do

   failures = 0
   do i = 1, n_cases
      if (differs(parallel_record(i), serial_record(i))) failures = failures + 1
   end do
   do i = 1, nwork
      call tmatrix_workspace_finalize(work(i))
   end do

   !! Regime census from the flags the serial pass recorded.  A harness that
   !! silently stopped covering a regime would still pass the bit comparison,
   !! so the coverage it claims is asserted rather than assumed.
   n_rayleigh  = count(serial_record%flag == 10)
   n_geometric = count(serial_record%flag == 20)
   n_tmatrix   = count(serial_record%flag == 0)
   if (n_rayleigh == 0)  error stop 'no case reached the Rayleigh dipole limit'
   if (n_tmatrix == 0)   error stop 'no case reached the T-matrix solver'
   if (n_geometric == 0) error stop 'no case reached the geometric optics limit'

   if (failures /= 0) error stop 'parallel spheroid_q or asymptotic-limit mismatch'
   write(*,'(a,i0,a,i0,a,i0,a,i0,a)') &
      'PASS: OpenMP spheroid_q and asymptotic limits (', nwork, ' workspaces; ', &
      n_rayleigh, ' Rayleigh, ', n_tmatrix, ' T-matrix, ', n_geometric, ' geometric optics)'

contains

   subroutine evaluate(work, icase, rec)
      !! One case: the regime-selecting entry point the host calls, followed by
      !! the two closed-form limits with their orientation-resolved outputs.
      type(tmatrix_workspace_t), intent(inout) :: work
      integer, intent(in) :: icase
      type(case_record_t), intent(out) :: rec
      type(tmatrix_options_t) :: options
      real(real64) :: a, lambda, nr, ki, aspect_ratio
      integer :: igeom, iabs, iasp

      igeom = mod(icase-1, n_geom) + 1
      iabs  = mod((icase-1)/n_geom, 2) + 1
      iasp  = (icase-1)/(2*n_geom) + 1
      a      = a_list(igeom)
      lambda = lam_list(igeom)
      ki     = absorption(iabs)
      aspect_ratio = aspect(iasp)
      if (mod(igeom+iabs, 2) == 0) then
         nr = 1.40_real64
      else
         nr = 1.70_real64
      end if

      options%tolerance    = tolerance
      options%ndgs         = ndgs
      options%shape        = spheroid_shape
      options%aspect_ratio = aspect_ratio

      call spheroid_q(work, a, lambda, nr, ki, options, &
                      rec%qext, rec%qabs, rec%qsca, rec%albedo, rec%asymmetry, &
                      rec%flag, rec%status)

      call rayleigh_limit(a, lambda, nr, ki, aspect_ratio, &
                          rec%ray_q(1), rec%ray_q(2), rec%ray_q(3), rec%ray_q(4), &
                          rec%ray_coeff(:,1), rec%ray_coeff(:,2), rec%ray_coeff(:,3), &
                          rec%ray_coeff(:,4), rec%ray_coeff(:,5), rec%ray_coeff(:,6), &
                          rec%ray_lmax, &
                          rec%ray_ori(:,1), rec%ray_ori(:,2), rec%ray_ori(:,3), rec%ray_ori(:,4))

      call geometric_optics_limit(a, lambda, nr, ki, aspect_ratio, &
                          rec%go_q(1), rec%go_q(2), rec%go_q(3), rec%go_q(4), &
                          rec%go_coeff(:,1), rec%go_coeff(:,2), rec%go_coeff(:,3), &
                          rec%go_coeff(:,4), rec%go_coeff(:,5), rec%go_coeff(:,6), &
                          rec%go_lmax, &
                          rec%go_ori(:,1), rec%go_ori(:,2), rec%go_ori(:,3), rec%go_ori(:,4))
   end subroutine evaluate

   logical function differs(actual, expected)
      type(case_record_t), intent(in) :: actual, expected
      differs = actual%flag   /= expected%flag   .or. &
                actual%status /= expected%status .or. &
                actual%ray_lmax /= expected%ray_lmax .or. &
                actual%go_lmax  /= expected%go_lmax  .or. &
                bits_differ([actual%qext, actual%qabs, actual%qsca, &
                             actual%albedo, actual%asymmetry], &
                            [expected%qext, expected%qabs, expected%qsca, &
                             expected%albedo, expected%asymmetry]) .or. &
                bits_differ(actual%ray_q, expected%ray_q) .or. &
                bits_differ(actual%go_q,  expected%go_q)  .or. &
                bits_differ(reshape(actual%ray_coeff, [n_coeff*6]), &
                            reshape(expected%ray_coeff, [n_coeff*6])) .or. &
                bits_differ(reshape(actual%go_coeff, [n_coeff*6]), &
                            reshape(expected%go_coeff, [n_coeff*6])) .or. &
                bits_differ(reshape(actual%ray_ori, [n_ori*4]), &
                            reshape(expected%ray_ori, [n_ori*4])) .or. &
                bits_differ(reshape(actual%go_ori, [n_ori*4]), &
                            reshape(expected%go_ori, [n_ori*4]))
   end function differs

   logical function bits_differ(actual, expected)
      !! Compare IEEE bit patterns.  The two runs execute the same code on the
      !! same inputs, so anything short of an identical pattern is state shared
      !! between threads, not a rounding difference.
      real(real64), intent(in) :: actual(:), expected(:)
      bits_differ = any(transfer(actual, 0_int64, size(actual)) /= &
                        transfer(expected, 0_int64, size(expected)))
   end function bits_differ

end program test_parallel_spheroid_optics
