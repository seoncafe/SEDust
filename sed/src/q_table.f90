module q_table_mod
   ! Reader for our run_tmatrix.x output
   !   tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat
   !
   ! File format (4 `#` header lines, then NA*NW data rows):
   !   lambda[um]  a_eff[um]  Q_ext  Q_abs  Q_sca  albedo  g  flag
   !
   ! Outer index in the loop is jw (lambda), inner index is ja (a_eff),
   ! per run_tmatrix.f90's loop nesting. Recommend running the FULL
   ! sweep (no 'test' arg) before invoking this module — the test-mode
   ! file has stride > 1 and breaks the NA*NW row count assumed here.
   !
   ! Also provides a simple log-linear interpolator in `a_eff` at fixed
   ! lambda index, which is the operation tau_check needs.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   implicit none
   private
   public :: load_q_table, interp_q_in_a
   public :: n_lam, n_aeff, lam_t, aeff_t, qext, qabs, qsca, albedo, gpar, flag

   integer, parameter :: wp = real64

   integer, parameter :: NA_DEF = 169   ! DH21 size grid length
   integer, parameter :: NW_DEF = 1129  ! DH21 wavelength grid length

   integer  :: n_aeff = 0, n_lam = 0
   real(wp), allocatable :: aeff_t(:), lam_t(:)            ! grid axes
   real(wp), allocatable :: qext(:,:), qabs(:,:), qsca(:,:)
   real(wp), allocatable :: albedo(:,:), gpar(:,:)
   integer,  allocatable :: flag(:,:)                       ! diagnostic
   real(wp), allocatable :: log_a(:)                        ! cached for interp

contains

   subroutine load_q_table(filename, na_in, nw_in)
      ! Reads a (NA, NW) Q table written by run_tmatrix.f90 in jw-outer,
      ! ja-inner order. Defaults na=169, nw=1129; pass other values when
      ! reading a different grid.
      character(len=*),  intent(in)           :: filename
      integer, optional, intent(in)           :: na_in, nw_in
      integer  :: u, ios, i, ja, jw, ifl
      real(wp) :: lam, ae, qe, qa, qs, w, g
      character(len=512) :: line

      n_aeff = NA_DEF;  if (present(na_in)) n_aeff = na_in
      n_lam  = NW_DEF;  if (present(nw_in)) n_lam  = nw_in

      if (allocated(aeff_t)) deallocate(aeff_t, lam_t, qext, qabs, qsca, albedo, gpar, flag, log_a)
      allocate(aeff_t(n_aeff), lam_t(n_lam))
      allocate(qext(n_lam, n_aeff), qabs(n_lam, n_aeff), qsca(n_lam, n_aeff))
      allocate(albedo(n_lam, n_aeff), gpar(n_lam, n_aeff), flag(n_lam, n_aeff))
      allocate(log_a(n_aeff))

      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(error_unit,'(a,a)') 'load_q_table: cannot open ', trim(filename)
         stop 1
      end if

      ! Skip leading `#` header lines.
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then
            write(error_unit,'(a)') 'load_q_table: unexpected EOF in header.'
            stop 1
         end if
         line = adjustl(line)
         if (len_trim(line) == 0)   cycle
         if (line(1:1) == '#')      cycle
         exit          ! line now holds first data record
      end do

      ! Process the first data line we already have.
      ja = 1; jw = 1
      read(line,*) lam, ae, qe, qa, qs, w, g, ifl
      lam_t(jw)         = lam
      aeff_t(ja)        = ae
      qext(jw, ja)      = qe
      qabs(jw, ja)      = qa
      qsca(jw, ja)      = qs
      albedo(jw, ja)    = w
      gpar(jw, ja)      = g
      flag(jw, ja)      = ifl

      do i = 2, n_aeff*n_lam
         read(u,*,iostat=ios) lam, ae, qe, qa, qs, w, g, ifl
         if (ios /= 0) then
            write(error_unit,'(a,i0,a,i0,a)') &
               'load_q_table: read error at row ', i, ' of ', n_aeff*n_lam, &
               '. Did you run the FULL sweep (./run_tmatrix.x) before tau_check?'
            stop 1
         end if
         ja = ja + 1
         if (ja > n_aeff) then
            ja = 1
            jw = jw + 1
            if (jw > n_lam) exit
         end if
         lam_t(jw)         = lam
         aeff_t(ja)        = ae
         qext(jw, ja)      = qe
         qabs(jw, ja)      = qa
         qsca(jw, ja)      = qs
         albedo(jw, ja)    = w
         gpar(jw, ja)      = g
         flag(jw, ja)      = ifl
      end do
      close(u)

      do ja = 1, n_aeff
         log_a(ja) = log10(aeff_t(ja))
      end do
   end subroutine load_q_table


   subroutine interp_q_in_a(jw, a_target, qe, qa, qs)
      ! Log-linear interpolation in a_eff at fixed lambda index jw.
      ! Clamps to grid edges (returns the boundary value if a_target
      ! is outside [aeff_t(1), aeff_t(n_aeff)]).
      integer,  intent(in)  :: jw
      real(wp), intent(in)  :: a_target
      real(wp), intent(out) :: qe, qa, qs
      integer  :: lo, hi, mid
      real(wp) :: x, t

      if (a_target <= aeff_t(1)) then
         qe = qext(jw, 1); qa = qabs(jw, 1); qs = qsca(jw, 1); return
      end if
      if (a_target >= aeff_t(n_aeff)) then
         qe = qext(jw, n_aeff); qa = qabs(jw, n_aeff); qs = qsca(jw, n_aeff); return
      end if
      x = log10(a_target)
      lo = 1; hi = n_aeff
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (log_a(mid) <= x) then
            lo = mid
         else
            hi = mid
         end if
      end do
      t  = (x - log_a(lo)) / (log_a(hi) - log_a(lo))
      qe = (1.0_wp - t) * qext(jw, lo) + t * qext(jw, hi)
      qa = (1.0_wp - t) * qabs(jw, lo) + t * qabs(jw, hi)
      qs = (1.0_wp - t) * qsca(jw, lo) + t * qsca(jw, hi)
   end subroutine interp_q_in_a

end module q_table_mod
