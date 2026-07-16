module pah_ld01_mod
   ! Reads Draine's tabulated carbonaceous (PAH/graphitic) Q_abs from the
   ! reference files PAHneu.31 / PAHion.31 (neutral / cation), which Draine
   ! and HD23 used. These are LD01 (Li & Draine 2001) for a < 0.05 um,
   ! lambda > 0.03 um (and DL84/LD93 elsewhere). Provides q_pah_ld01(icharge,
   ! a_um, lam_um) -> Q_abs by bilinear interpolation in (log a, log lambda),
   ! so the SED solver can use the *reference* cross sections directly instead
   ! of our DL07 qpah implementation. Set use_ld01_pah_xsec=.true. to enable.
   use constants, only: wp
   implicit none
   private
   public :: q_pah_ld01, load_pah_ld01, use_ld01_pah_xsec

   logical, save :: use_ld01_pah_xsec = .false.   ! toggled by the driver

   integer, parameter :: NRAD = 31, NWAV = 1201
   character(len=*), parameter :: F_NEU = '../data/dielectric/PAHneu.31'
   character(len=*), parameter :: F_ION = '../data/dielectric/PAHion.31'

   logical,  save :: loaded = .false.
   real(wp), save :: la(NRAD)              ! log(radius/um), ascending
   real(wp), save :: lw(NWAV)              ! log(lambda/um), ascending
   real(wp), save :: lqn(NWAV,NRAD), lqi(NWAV,NRAD)   ! log(Q_abs) neu/ion

contains

   subroutine load_pah_ld01(ok)
      ! Optional ok: absent -> stop on error as before; present -> return
      ! .false. (leaving loaded=.false.) instead of stopping. When ok is
      ! present it is forwarded to read_one so the reader stays silent and
      ! reports through the flag; when absent read_one keeps its own message.
      logical, optional, intent(out) :: ok
      logical :: ok1
      if (present(ok)) ok = .true.
      if (loaded) return
      if (present(ok)) then
         call read_one(F_NEU, lqn, .true.,  ok1)
         if (.not. ok1) then;  ok = .false.;  return;  end if
         call read_one(F_ION, lqi, .false., ok1)
         if (.not. ok1) then;  ok = .false.;  return;  end if
      else
         call read_one(F_NEU, lqn, .true.)
         call read_one(F_ION, lqi, .false.)
      end if
      loaded = .true.
   end subroutine load_pah_ld01

   subroutine read_one(fname, lq, set_grids, ok)
      character(len=*), intent(in)  :: fname
      real(wp),         intent(out) :: lq(NWAV,NRAD)
      logical,          intent(in)  :: set_grids
      logical, optional, intent(out) :: ok
      integer  :: u, ios, ir, iw, ist
      real(wp) :: ww, qe, qa
      character(len=256) :: line
      real(wp) :: araw(NRAD), wraw(NWAV), qraw(NWAV,NRAD)

      if (present(ok)) ok = .true.
      open(newunit=u, file=fname, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(ok)) then
            ok = .false.;  return
         else
            write(*,'(a,a)') ' pah_ld01: cannot open ', trim(fname); stop 1
         end if
      end if
      ir = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (index(line,'radius(micron)') > 0) then
            ir = ir + 1
            read(line,*) araw(ir)
            read(u,'(a)') line                 ! skip the column header
            do iw = 1, NWAV
               read(u,'(a)') line
               ! parse only w (item 1), Q_ext (item 2), Q_abs (item 3);
               ! Q_sca/g (items 4-5) can have malformed exponents -> skip them.
               read(line,*,iostat=ist) ww, qe, qa
               if (ist /= 0) then
                  if (present(ok)) then
                     close(u);  ok = .false.;  return
                  else
                     write(*,'(a,i0,a,i0,a)') ' pah_ld01: parse fail ir=', ir, &
                          ' iw=', iw, ' line=['//trim(line)//']'; stop 1
                  end if
               end if
               qraw(iw,ir) = qa
               if (ir == 1) wraw(iw) = ww
            end do
         end if
      end do
      close(u)

      ! wraw is descending (1000 -> 1e-3 um); reverse to ascending and
      ! store log values. araw is ascending (3.16e-4 -> 1e-2 um).
      do iw = 1, NWAV
         do ir = 1, NRAD
            lq(iw,ir) = log(max(qraw(NWAV-iw+1,ir), tiny(0.0_wp)))
         end do
      end do
      if (set_grids) then
         do ir = 1, NRAD
            la(ir) = log(araw(ir))
         end do
         do iw = 1, NWAV
            lw(iw) = log(wraw(NWAV-iw+1))
         end do
      end if
   end subroutine read_one

   function q_pah_ld01(icharge, a_um, lam_um) result(Qabs)
      integer,  intent(in) :: icharge        ! 0 = neutral, 1 = cation
      real(wp), intent(in) :: a_um, lam_um
      real(wp) :: Qabs
      real(wp) :: xa, xw, ta, tw, q00, q01, q10, q11
      integer  :: ia, iw

      if (.not. loaded) call load_pah_ld01()

      ! clamp to the table range (a>100AA contributes <0.1% of PAH area)
      xa = min(max(log(a_um),   la(1)), la(NRAD))
      xw = min(max(log(lam_um), lw(1)), lw(NWAV))
      ia = bracket(la, NRAD, xa)
      iw = bracket(lw, NWAV, xw)
      ta = (xa - la(ia)) / (la(ia+1) - la(ia))
      tw = (xw - lw(iw)) / (lw(iw+1) - lw(iw))

      if (icharge == 0) then
         q00 = lqn(iw,ia);     q10 = lqn(iw,ia+1)
         q01 = lqn(iw+1,ia);   q11 = lqn(iw+1,ia+1)
      else
         q00 = lqi(iw,ia);     q10 = lqi(iw,ia+1)
         q01 = lqi(iw+1,ia);   q11 = lqi(iw+1,ia+1)
      end if
      Qabs = exp( (1-tw)*((1-ta)*q00 + ta*q10) + tw*((1-ta)*q01 + ta*q11) )
   end function q_pah_ld01

   pure function bracket(x, n, xv) result(i)
      integer,  intent(in) :: n
      real(wp), intent(in) :: x(n), xv
      integer :: i, lo, hi, mid
      lo = 1; hi = n
      do while (hi - lo > 1)
         mid = (lo + hi)/2
         if (x(mid) <= xv) then; lo = mid; else; hi = mid; end if
      end do
      i = lo
   end function bracket

end module pah_ld01_mod
