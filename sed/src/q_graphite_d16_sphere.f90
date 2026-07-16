module q_graphite_d16_sphere_mod
   ! D16 turbostratic graphite (Draine 2016, MG EMT) on *spheres*.
   ! Reads the content of data/dielectric/callqcomp_D16MGemt.gz, which we
   ! keep decompressed as data/dielectric/q_D16graphite.dat.
   !
   ! Used to isolate the sphere-vs-spheroid contribution in the
   ! qpah.f90 xi-blend: comparing this against q_graphite_d16_mod
   ! (b/a = 1.4 oblate spheroid, qlib_gra_D16MGemt_1.400) at fixed
   ! material (D16 MG EMT) reveals how much of the D03-sphere ->
   ! D16-spheroid sub-mm jump is shape vs material.
   !
   ! File layout (q_D16graphite.dat):
   !   line 1            : header '41 3501   # NRAD NWAV (...)'
   !   line 2            : 41 radii (um), ascending 1e-3 to 1e+1
   !   line 3            : 3501 wavelengths (um), DESCENDING 1e+4 to 1e-3
   !   lines 4-44        : 41 rows of 3501 Q_abs values; row i = radius i

   use constants, only: wp
   implicit none
   private
   public :: q_graphite_d16_sphere_abs, load_q_graphite_d16_sphere

   character(len=*), parameter :: F_D16S = '../data/dielectric/q_D16graphite.dat'
   integer,  parameter :: NA_S = 41
   integer,  parameter :: NW_S = 3501

   logical  :: loaded = .false.
   real(wp) :: a_grid(NA_S),  log_a_grid(NA_S)
   real(wp) :: w_grid(NW_S),  log_w_grid(NW_S)        ! descending in file
   real(wp) :: Q_tab(NA_S, NW_S)                      ! Q_abs(jrad, jwave)

contains

   subroutine load_q_graphite_d16_sphere(ok)
      ! Optional ok: absent -> stop on error as before; present -> return
      ! .false. (leaving loaded=.false.) instead of stopping.
      logical, optional, intent(out) :: ok
      integer  :: u, ios, ja
      character(len=512) :: hdr

      if (present(ok)) ok = .true.
      open(newunit=u, file=F_D16S, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(ok)) then
            ok = .false.;  return
         else
            write(*,'(a,a)') 'q_graphite_d16_sphere: cannot open ', F_D16S
            stop 1
         end if
      end if
      read(u, '(a)') hdr             ! '41 3501   # NRAD NWAV ...'
      read(u, *)     a_grid          ! 41 radii (um)
      read(u, *)     w_grid          ! 3501 wavelengths (um) -- descending
      do ja = 1, NA_S
         read(u, *) Q_tab(ja, :)
      end do
      close(u)

      log_a_grid = log(a_grid)
      log_w_grid = log(w_grid)
      loaded = .true.
   end subroutine load_q_graphite_d16_sphere


   subroutine q_graphite_d16_sphere_abs(agrain, lambda, Qabs)
      ! Bilinear interp in (log a, log lambda); wavelengths descending
      ! so we flip the bracket sense.
      real(wp), intent(in)  :: agrain      ! [um]
      real(wp), intent(in)  :: lambda      ! [um]
      real(wp), intent(out) :: Qabs
      real(wp) :: la, lw, ta, tw, q00, q01, q10, q11
      integer  :: ia, iw

      if (.not. loaded) call load_q_graphite_d16_sphere()

      la = log(agrain)
      lw = log(lambda)
      call bracket_asc (log_a_grid, NA_S, la, ia, ta)
      call bracket_desc(log_w_grid, NW_S, lw, iw, tw)

      q00 = Q_tab(ia,   iw  )
      q10 = Q_tab(ia+1, iw  )
      q01 = Q_tab(ia,   iw+1)
      q11 = Q_tab(ia+1, iw+1)
      Qabs = (1.0_wp-ta)*(1.0_wp-tw)*q00 + ta*(1.0_wp-tw)*q10 &
           + (1.0_wp-ta)*tw       *q01 + ta*tw       *q11
   end subroutine q_graphite_d16_sphere_abs


   subroutine bracket_asc(grid, n, x, i, t)
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: grid(n), x
      integer,  intent(out) :: i
      real(wp), intent(out) :: t
      integer :: lo, hi, mid

      if (x <= grid(1)) then;  i = 1;   t = 0.0_wp;  return; end if
      if (x >= grid(n)) then;  i = n-1; t = 1.0_wp;  return; end if
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (grid(mid) <= x) then; lo = mid; else; hi = mid; end if
      end do
      i = lo
      t = (x - grid(i)) / (grid(i+1) - grid(i))
   end subroutine bracket_asc


   subroutine bracket_desc(grid, n, x, i, t)
      ! grid descending: grid(1) >= grid(2) >= ... >= grid(n).
      ! Return i such that grid(i) >= x >= grid(i+1), and
      ! t = (grid(i) - x) / (grid(i) - grid(i+1)) so that
      ! value(x) = (1-t) value(i) + t value(i+1).
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: grid(n), x
      integer,  intent(out) :: i
      real(wp), intent(out) :: t
      integer :: lo, hi, mid

      if (x >= grid(1)) then;  i = 1;   t = 0.0_wp;  return; end if
      if (x <= grid(n)) then;  i = n-1; t = 1.0_wp;  return; end if
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (grid(mid) >= x) then; lo = mid; else; hi = mid; end if
      end do
      i = lo
      t = (grid(i) - x) / (grid(i) - grid(i+1))
   end subroutine bracket_desc

end module q_graphite_d16_sphere_mod
