module q_graphite_d16_mod
   ! D16 turbostratic graphite (Draine 2016, MG EMT) on oblate spheroids
   ! with b/a = 1.4. Random-orientation-averaged Q_abs(a, lambda), read
   ! from Draine's precomputed qlib_gra_D16MGemt_1.400 (T-matrix output).
   !
   ! Drop-in replacement for q_graphite_mod::q_graphite_abs in qpah.f90's
   ! HD23 eq. 15-16 / DL07 eq. 5-7 PAH-to-graphite transition.
   !
   ! File layout (qlib_gra_D16MGemt_1.400, decompressed):
   !   line 1            : title 'D16 graphite: turbostratic, MG_1 EMT'
   !   line 2            : NRAD NWAV (= 168 1008, meaning indices 0..NRAD)
   !   line 3            : AXRAT (b/a = 1.4)
   !   line 4            : 'rad(0)...rad(NRAD) follow:'
   !   subsequent lines  : (NRAD+1) = 169 radii in um, 8 per line
   !   marker line       : 'wave(0)...wav(NWAV) follow:'
   !   subsequent lines  : (NWAV+1) = 1009 wavelengths in um, 8 per line
   !   marker line       : 'Qabs(jori=1-3,jrad=0-NRAD,jwave=0-NWAV) follow:'
   !   data stream       : 3*(NRAD+1)*(NWAV+1) values, jori innermost,
   !                       jrad middle, jwave outermost; 8 per line until
   !                       the final line which holds the remainder.
   !   then              : Qext block, Qsca block (we ignore them).
   !
   ! jori convention (Mishchenko):
   !   1 = random-orientation average  <- what ksi-blend wants
   !   2 = k perp a, E parallel to a   (E along symmetry axis)
   !   3 = k perp a, E perp a

   use constants, only: wp
   implicit none
   private
   public :: q_graphite_d16_abs, load_q_graphite_d16

   character(len=*), parameter :: F_D16 = '../data/dielectric/qlib_gra_D16MGemt_1.400'
   integer,  parameter :: NA_D = 169
   integer,  parameter :: NW_D = 1009

   logical  :: loaded = .false.
   real(wp) :: a_grid(NA_D),      log_a_grid(NA_D)
   real(wp) :: w_grid(NW_D),      log_w_grid(NW_D)
   real(wp) :: Q_rand(NA_D, NW_D)            ! Q_abs(jori=1)

contains

   subroutine load_q_graphite_d16(ok)
      ! Optional ok: absent -> stop on error as before; present -> return
      ! .false. (leaving loaded=.false.) instead of stopping.
      logical, optional, intent(out) :: ok
      integer  :: u, ios, nrad_h, nwav_h
      character(len=256) :: hdr
      real(wp), allocatable :: Q_all(:,:,:)
      real(wp) :: axrat

      if (present(ok)) ok = .true.
      open(newunit=u, file=F_D16, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(ok)) then
            ok = .false.;  return
         else
            write(*,'(a,a)') 'q_graphite_d16: cannot open ', F_D16
            stop 1
         end if
      end if

      read(u, '(a)') hdr                    ! title
      read(u, *)     nrad_h, nwav_h         ! 168 1008 (indices 0..N)
      if (nrad_h+1 /= NA_D .or. nwav_h+1 /= NW_D) then
         if (present(ok)) then
            close(u);  ok = .false.;  return
         else
            write(*,'(a,2i6)') 'q_graphite_d16: dim mismatch ', nrad_h, nwav_h
            stop 1
         end if
      end if
      read(u, *)     axrat                  ! 1.4
      read(u, '(a)') hdr                    ! 'rad(0)...rad(NRAD) follow:'
      read(u, *)     a_grid                 ! NA_D radii (um)
      read(u, '(a)') hdr                    ! 'wave(0)...wav(NWAV) follow:'
      read(u, *)     w_grid                 ! NW_D wavelengths (um)
      read(u, '(a)') hdr                    ! 'Qabs(...) follow:'

      ! Stream-read Qabs(jori=1..3, jrad=0..NRAD, jwave=0..NWAV).
      ! Fortran column-major layout has jori innermost -> matches file order.
      allocate(Q_all(3, NA_D, NW_D))
      read(u, *) Q_all
      Q_rand = Q_all(1, :, :)
      deallocate(Q_all)
      close(u)

      log_a_grid = log(a_grid)
      log_w_grid = log(w_grid)
      loaded = .true.
   end subroutine load_q_graphite_d16


   subroutine q_graphite_d16_abs(agrain, lambda, Qabs)
      ! Bilinear interpolation in (log a, log lambda) of jori=1 Q_abs.
      ! Out-of-grid inputs are clamped to the nearest grid edge.
      real(wp), intent(in)  :: agrain      ! [um]
      real(wp), intent(in)  :: lambda      ! [um]
      real(wp), intent(out) :: Qabs        ! dimensionless
      real(wp) :: la, lw, ta, tw, q00, q01, q10, q11
      integer  :: ia, iw

      if (.not. loaded) call load_q_graphite_d16()

      la = log(agrain)
      lw = log(lambda)
      call bracket_asc(log_a_grid, NA_D, la, ia, ta)
      call bracket_asc(log_w_grid, NW_D, lw, iw, tw)

      q00 = Q_rand(ia,   iw  )
      q10 = Q_rand(ia+1, iw  )
      q01 = Q_rand(ia,   iw+1)
      q11 = Q_rand(ia+1, iw+1)
      Qabs = (1.0_wp-ta)*(1.0_wp-tw)*q00 + ta*(1.0_wp-tw)*q10 &
           + (1.0_wp-ta)*tw       *q01 + ta*tw       *q11
   end subroutine q_graphite_d16_abs


   subroutine bracket_asc(grid, n, x, i, t)
      ! For an ascending grid, return i such that grid(i) <= x <= grid(i+1)
      ! and t in [0,1] = (x - grid(i)) / (grid(i+1) - grid(i)). Clamps
      ! at edges (i = 1, t = 0 below; i = n-1, t = 1 above).
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: grid(n), x
      integer,  intent(out) :: i
      real(wp), intent(out) :: t
      integer :: lo, hi, mid

      if (x <= grid(1)) then
         i = 1;   t = 0.0_wp;  return
      end if
      if (x >= grid(n)) then
         i = n-1; t = 1.0_wp;  return
      end if
      lo = 1
      hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (grid(mid) <= x) then
            lo = mid
         else
            hi = mid
         end if
      end do
      i = lo
      t = (x - grid(i)) / (grid(i+1) - grid(i))
   end subroutine bracket_asc

end module q_graphite_d16_mod
