module pah_ioniz_mod
   ! ----------------------------------------------------------------------
   ! Diffuse-ISM PAH ionization fraction, ported from B.T. Draine's
   ! Fortran-77 WD01b grain-charging stack (irem/src).
   !
   ! Public entry point:
   !
   !    f_ion = pah_ionfrac(a_um, U)
   !
   ! returns the fraction of PAHs of (graphite-equivalent) radius a_um
   ! [micron] that are NOT neutral, in the average diffuse ISM
   ! (CNM 43% + WNM 43% + WIM 14%), with the ISRF scaled to U*MMP83.
   ! This reproduces PAH_CHARGING_DISM(A,U,FRAC_IONIZATION).
   !
   ! This module follows Draine's grain-charging method (WD01b). Its
   ! procedures correspond to that method's steps as follows:
   !    PAH_CHARGING_DISM  -> pah_ionfrac (this module)
   !    CHARGE             -> charge (subroutine)
   !    DJPEDE, DGPEDE     -> djpede, dgpede (integrands)
   !    THRESHOLDS         -> thresholds
   !    PEYIELD            -> peyield
   !    DQG32              -> dqg32 (32-pt Gauss)
   !    RADFLD (MODE=4)    -> radfld_mmp
   !
   ! The Fortran-77 PEYIELD calls INDEX(...) only to obtain the graphite
   ! imaginary refractive index Im(n) for the small-particle yield
   ! enhancement factor Y1.  The QCOMP cross-section package (2798 lines)
   ! is NOT ported; instead the photon-weighting cross section Q_abs and
   ! the Y1 Im(n) are both computed from the Draine 2003 graphite
   ! dielectric tables (1/3 parallel + 2/3 perpendicular, with the
   ! size-dependent free-electron Drude term), exactly as in the project's
   ! q_graphite.f90.  See the "Deviations" note at the bottom of this file.
   !
   ! COMMON/JPECOM/ state is carried in module-level SAVEd variables.
   ! DQG32's EXTERNAL integrand is passed through an abstract interface.
   !
   ! Self-contained: depends only on `constants` (for wp).  All numerical
   ! helpers (Mie, linear interpolation, Planck) are private to this module.
   ! ----------------------------------------------------------------------
   use constants, only: wp
   implicit none
   private
   public :: pah_ionfrac

   ! --- numerical constants ------------------------------------------------
   real(wp), parameter :: PI = 3.141592653589793238462643383279502884197_wp

   ! --- graphite dielectric tables (D03) -----------------------------------
   character(len=*), parameter :: F_CPA = '../data/dielectric/index_CpaD03'
   character(len=*), parameter :: F_CPE = '../data/dielectric/index_CpeD03'
   integer,  parameter :: NCPA = 387
   integer,  parameter :: NCPE = 384
   real(wp), parameter :: T_GRAIN = 20.0_wp   ! K, Draine default (TG_COM=20)

   ! PAH<->graphite transition constants (QCOMP ICOMP=11-16)
   real(wp), parameter :: A_T   = 0.0050_wp   ! transition radius [um]
   real(wp), parameter :: FGMIN = 0.01_wp     ! minimum graphitic contribution

   logical  :: initialized = .false.
   ! base tables (bound-electron eps), wavelengths in micron (descending)
   real(wp) :: cpa_eV(NCPA), cpa_eps1(NCPA), cpa_eps2(NCPA), cpa_wavl(NCPA)
   real(wp) :: cpe_eV(NCPE), cpe_eps1(NCPE), cpe_eps2(NCPE), cpe_wavl(NCPE)
   ! cached (n,k) including free-electron term, rebuilt when radius changes
   real(wp) :: cached_a = -1.0_wp
   real(wp) :: cpa_n(NCPA), cpa_k(NCPA)
   real(wp) :: cpe_n(NCPE), cpe_k(NCPE)

   ! --- "COMMON/JPECOM/" state passed into the DQG32 integrands ------------
   integer  :: rfion_com, rftyp_com, icomp_com, jz_com
   real(wp) :: rfchi_com, rfpar_com, rfr_v_com, rftau_com
   real(wp) :: ar_com, tg_com, ip_com, emin_com, epdt_com, epet_com
   real(wp) :: wbulk_com

   ! abstract interface for the integrand passed to dqg32
   abstract interface
      function integrand(e) result(f)
         import :: wp
         real(wp), intent(in) :: e
         real(wp) :: f
      end function integrand
   end interface

contains

   ! ======================================================================
   !  Public driver: PAH_CHARGING_DISM
   ! ======================================================================
   function pah_ionfrac(a_um, U) result(f_ion)
      ! a_um : graphite-equivalent radius [micron]
      ! U    : ISRF intensity relative to MMP83
      real(wp), intent(in) :: a_um, U
      real(wp) :: f_ion

      integer, parameter :: IZMX = 55
      integer :: icomp, rftyp, rfion
      real(wp) :: a, rfchi, rfpar, rfr_v, rftau
      real(wp) :: frac_cnm, frac_wnm, frac_wim
      real(wp) :: fz(-IZMX:IZMX)
      integer  :: izmin, izmax, izmode

      a = a_um * 1.0e-4_wp   ! radius in cm

      icomp = 11             ! PAH/graphitic material
      rftyp = 4              ! MMP-type radiation field
      rfpar = 0.0_wp         ! not used when rftyp=4
      rfchi = 1.23_wp * U    ! intensity at 1000A rel. to Habing
      rfion = 0              ! cutoff at Lyman continuum
      rfr_v = 3.1_wp         ! R_V = A_V/E(B-V)
      rftau = 0.0_wp         ! tau_V

      if (.not. initialized) call init_module()

      ! --- CNM: T=100, nH=30, n(H+)=.0012nH, n(M+)=.0003nH, M+=12 mH -------
      call charge(IZMX, a, icomp, rftyp, rfion, rfchi, rfpar, rfr_v, rftau, &
                  100.0_wp, 30.0_wp, 0.0012_wp, 0.0003_wp, 12.0_wp, &
                  fz, izmin, izmax, izmode)
      frac_cnm = sum_nonneutral(fz, izmin, izmax, IZMX)

      ! --- WNM: T=6000, nH=0.4, n(H+)=0.1nH, n(M+)=.0003nH, M+=12 mH -------
      call charge(IZMX, a, icomp, rftyp, rfion, rfchi, rfpar, rfr_v, rftau, &
                  6000.0_wp, 0.4_wp, 0.1_wp, 0.0003_wp, 12.0_wp, &
                  fz, izmin, izmax, izmode)
      frac_wnm = sum_nonneutral(fz, izmin, izmax, IZMX)

      ! --- WIM: T=8000, nH=0.1, n(H+)=0.99nH, n(M+)=.001nH, M+=12 mH -------
      call charge(IZMX, a, icomp, rftyp, rfion, rfchi, rfpar, rfr_v, rftau, &
                  8000.0_wp, 0.1_wp, 0.99_wp, 0.001_wp, 12.0_wp, &
                  fz, izmin, izmax, izmode)
      frac_wim = sum_nonneutral(fz, izmin, izmax, IZMX)

      f_ion = 0.43_wp*frac_cnm + 0.43_wp*frac_wnm + 0.14_wp*frac_wim
   end function pah_ionfrac


   function sum_nonneutral(fz, izmin, izmax, izmx) result(s)
      integer,  intent(in) :: izmin, izmax, izmx
      real(wp), intent(in) :: fz(-izmx:izmx)
      real(wp) :: s
      integer  :: iz
      s = 0.0_wp
      do iz = izmin, izmax
         if (iz /= 0) s = s + fz(iz)
      end do
   end function sum_nonneutral


   ! ======================================================================
   !  CHARGE : Weingartner & Draine (1999) charge-distribution solver
   !  Note: x_h = n(H+)/n_H ; x_m = n(M+)/n_H ; m_i = effective metal mass.
   !        Here x_h and x_m enter as the (H+,M+)/n_H ratios passed in.
   ! ======================================================================
   subroutine charge(jzmx, a, icomp, rftyp, rfion, rfchi, rfpar, &
                     rfr_v, rftau, t, n_h, x_h, x_m, m_i, &
                     fz, izmin, izmax, izmode)
      integer,  intent(in)  :: jzmx, icomp, rftyp, rfion
      real(wp), intent(in)  :: a, rfchi, rfpar, rfr_v, rftau
      real(wp), intent(in)  :: t, n_h, x_h, x_m, m_i
      real(wp), intent(out) :: fz(-jzmx:jzmx)
      integer,  intent(out) :: izmin, izmax, izmode

      real(wp), parameter :: ARGMAX = 50.0_wp, FMIN = 1.0e-10_wp
      integer  :: jjz, jz, jzmin, jzmax
      real(wp) :: aphys, djpe, ephmax, ephmx, fac
      real(wp) :: j0, je0, ji0, n_m
      real(wp) :: se0, sep, ssum, tau, term, wbulk
      real(wp) :: ip(-jzmx:jzmx), emin(-jzmx:jzmx)
      real(wp) :: epdt(-jzmx:jzmx), epet(-jzmx:jzmx)
      real(wp) :: je(-jzmx:jzmx), ji(-jzmx:jzmx), jpe(-jzmx:jzmx)
      real(wp) :: theta(1:jzmx)

      ! initialize arrays
      fz   = 0.0_wp
      ip   = 0.0_wp; emin = 0.0_wp; epdt = 0.0_wp; epet = 0.0_wp
      je   = 0.0_wp; ji   = 0.0_wp; jpe  = 0.0_wp
      izmin = 0; izmax = 0; izmode = 0

      tg_com = T_GRAIN

      do jz = 1, jzmx
         theta(jz) = jz / (1.0_wp + 1.0_wp/sqrt(real(jz, wp)))
      end do

      ! work function (graphite); aphys = a (no APHYS/A distinction)
      wbulk = 4.4_wp
      aphys = a
      wbulk_com = wbulk

      ! tau = a k T / e^2  (Draine & Sutin 1987)
      tau = aphys*1.38e-16_wp*t/(4.803e-10_wp)**2

      ! most negative charge state (autoionization limit)
      jzmin = -jzmx
      do jz = -1, -jzmx, -1
         call thresholds(icomp, jz, a, ip(jz), emin(jz), epdt(jz), epet(jz))
         if ((ip(jz) + emin(jz)) < 0.0_wp) then
            jzmin = jz + 1
            exit
         end if
      end do

      ! maximum photon energy (eV)
      ephmax = 13.6_wp
      if (rftyp == 0 .and. rfion == 1) then
         ephmax = 100.0_wp
      else if (rftyp == 3 .and. rfion == 1) then
         ephmax = 10.0_wp*8.6171e-5_wp*rfpar
      end if

      ! pass radiation-field/grain state to integrands
      icomp_com = icomp
      rfion_com = rfion
      rftyp_com = rftyp
      rfchi_com = rfchi
      rfpar_com = rfpar
      rfr_v_com = rfr_v
      rftau_com = rftau
      ar_com    = a

      ! collision rates for neutral grain (no induced dipole)
      j0  = 1.9517e6_wp*n_h*sqrt(t)*aphys**2
      je0 = j0*(x_h + x_m)
      ji0 = j0*(x_h + x_m/sqrt(m_i))/sqrt(1836.1_wp)

      ! number of carbon atoms (graphitic)
      n_m = 468.0_wp*(1.0e7_wp*a)**3

      sep = 0.5_wp*(1.0_wp - exp(-1.0e7_wp*aphys))
      se0 = sep/(1.0_wp + exp(20.0_wp - n_m))

      if (jzmin < 0) then
         ! consider charge states -1 and 0 to find peak direction
         jz_com   = -1
         ip_com   = ip(-1)
         emin_com = emin(-1)
         epdt_com = epdt(-1)
         epet_com = epet(-1)

         if (ephmax > 10.0_wp*(1.0_wp + epdt_com)) then
            ephmx = sqrt(ephmax*(1.0_wp + epdt_com))
            call dqg32(epdt_com, ephmx,  djpede, jpe(-1))
            call dqg32(ephmx,  ephmax,   djpede, djpe)
            jpe(-1) = jpe(-1) + djpe
         else
            if (ephmax > epdt_com) then
               call dqg32(epdt_com, ephmax, djpede, jpe(-1))
            else
               jpe(-1) = 0.0_wp
            end if
         end if

         ji(-1) = ji0*(1.0_wp + 1.0_wp/tau)*(1.0_wp + sqrt(2.0_wp/(tau + 2.0_wp)))

         je(0) = je0*se0*(1.0_wp + sqrt(0.5_wp*PI/tau))
         if (ip(-1) < 0.0_wp) then
            if (ip(-1)/t > -1.0e-2_wp) then
               je(0) = je(0)*exp(1.1604e4_wp*ip(-1)/t)
            else
               je(0) = 0.0_wp
            end if
         end if

         fz(0)  = 1.0_wp
         fz(-1) = fz(0)*je(0)/(jpe(-1) + ji(-1))
         ssum   = fz(0) + fz(-1)
      else
         fz(0)  = 1.0_wp
         fz(-1) = 0.0_wp
         ssum   = 1.0_wp
      end if

      if (fz(0) > fz(-1)) then
         ! ---------------- peak at positive Z -----------------------------
         izmax = 0
         do jz = 0, jzmx-1
            call thresholds(icomp, jz, a, ip(jz), emin(jz), epdt(jz), epet(jz))
            jz_com   = jz
            ip_com   = ip(jz)
            emin_com = emin(jz)
            epdt_com = epdt(jz)
            epet_com = epet(jz)
            if (ephmax > 10.0_wp*(1.0_wp + epet_com)) then
               ephmx = sqrt(ephmax*(1.0_wp + epet_com))
               call dqg32(epet_com, ephmx,  djpede, jpe(jz))
               call dqg32(ephmx,   ephmax,  djpede, djpe)
               jpe(jz) = jpe(jz) + djpe
            else
               if (ephmax > epet_com) then
                  call dqg32(epet_com, ephmax, djpede, jpe(jz))
               else
                  jpe(jz) = 0.0_wp
               end if
            end if

            ! electrons attracted
            je(jz+1) = je0*sep*(1.0_wp + real(jz+1, wp)/tau)* &
                       (1.0_wp + sqrt(2.0_wp/(tau + real(2*jz+2, wp))))
            if (jz == 0) then
               ji(0) = ji0*(1.0_wp + sqrt(0.5_wp*PI/tau))
            else
               term = theta(jz)/tau
               if (term < ARGMAX) then
                  ji(jz) = ji0*(1.0_wp + 1.0_wp/sqrt(4.0_wp*tau + real(3*jz-1, wp)))**2 &
                           *exp(-term)
               else
                  ji(jz) = 0.0_wp
               end if
            end if
            fz(jz+1) = fz(jz)*(jpe(jz) + ji(jz))/je(jz+1)

            if (fz(jz+1) > 1.0e50_wp) then
               fac = 1.0_wp/fz(jz+1)
               do jjz = -1, jz+1
                  fz(jjz) = fac*fz(jjz)
               end do
               ssum = fac*ssum
            end if
            if (fz(jz+1) > 0.0_wp) then
               izmax = jz+1
               ssum  = ssum + fz(jz+1)
            end if
            if (fz(jz+1) < FMIN*ssum) exit
         end do

         izmin = -1
         if (jzmin < -1) then
            do jz = -2, jzmin, -1
               jz_com   = jz
               emin_com = emin(jz)
               epdt_com = epdt(jz)
               epet_com = epet(jz)
               if (ephmax > 10.0_wp*(epdt(jz) + 1.0_wp)) then
                  ephmx = sqrt((1.0_wp + epdt(jz))*ephmax)
                  call dqg32(epdt_com, ephmx,  djpede, jpe(jz))
                  call dqg32(ephmx,   ephmax,  djpede, djpe)
                  jpe(jz) = jpe(jz) + djpe
               else
                  if (ephmax > epdt_com) then
                     call dqg32(epdt_com, ephmax, djpede, jpe(jz))
                  else
                     jpe(jz) = 0.0_wp
                  end if
               end if

               ji(jz) = ji0*(1.0_wp - real(jz, wp)/tau)* &
                        (1.0_wp + sqrt(2.0_wp/(tau - real(2*jz, wp))))

               term = theta(-jz-1)/tau
               if (term > ARGMAX) then
                  izmin = jz+1
                  exit
               end if
               je(jz+1) = je0*se0* &
                          (1.0_wp + 1.0_wp/sqrt(4.0_wp*tau - real(3*jz+3, wp)))**2 &
                          *exp(-term)
               fz(jz) = fz(jz+1)*je(jz+1)/(jpe(jz) + ji(jz))
               izmin = jz
               ssum = ssum + fz(jz)
               if (fz(jz) < FMIN*ssum) exit
            end do
         end if
      else
         ! ---------------- peak at negative Z -----------------------------
         izmin = -1
         if (jzmin < -1) then
            do jz = -2, jzmin, -1
               jz_com   = jz
               emin_com = emin(jz)
               epdt_com = epdt(jz)
               epet_com = epet(jz)
               if (ephmax > 10.0_wp*(epdt(jz) + 1.0_wp)) then
                  ephmx = sqrt((1.0_wp + epdt(jz))*ephmax)
                  call dqg32(epdt_com, ephmx,  djpede, jpe(jz))
                  call dqg32(ephmx,   ephmax,  djpede, djpe)
                  jpe(jz) = jpe(jz) + djpe
               else
                  if (ephmax > epdt_com) then
                     call dqg32(epdt_com, ephmax, djpede, jpe(jz))
                  else
                     jpe(jz) = 0.0_wp
                  end if
               end if

               ji(jz) = ji0*(1.0_wp - real(jz, wp)/tau)* &
                        (1.0_wp + sqrt(2.0_wp/(tau - real(2*jz, wp))))

               term = theta(-jz-1)/tau
               if (term > ARGMAX) then
                  izmin = jz+1
                  exit
               end if
               je(jz+1) = je0*se0* &
                          (1.0_wp + 1.0_wp/sqrt(4.0_wp*tau - real(3*jz+3, wp)))**2 &
                          *exp(-term)
               fz(jz) = fz(jz+1)*je(jz+1)/(jpe(jz) + ji(jz))
               izmin = jz
               ssum = ssum + fz(jz)

               if (fz(jz) > 1.0e50_wp) then
                  fac = 1.0_wp/fz(jz)
                  do jjz = jz, 0
                     fz(jjz) = fac*fz(jjz)
                  end do
                  ssum = fac*ssum
               end if
               if (fz(jz) < FMIN*ssum) exit
            end do
         end if

         izmax = 0
         do jz = 0, jzmx-1
            call thresholds(icomp, jz, a, ip(jz), emin(jz), epdt(jz), epet(jz))
            jz_com   = jz
            emin_com = emin(jz)
            epdt_com = epdt(jz)
            epet_com = epet(jz)
            if (ephmax > 10.0_wp*(1.0_wp + epet_com)) then
               ephmx = sqrt(ephmax*(1.0_wp + epet_com))
               call dqg32(epet_com, ephmx,  djpede, jpe(jz))
               call dqg32(ephmx,   ephmax,  djpede, djpe)
               jpe(jz) = jpe(jz) + djpe
            else
               if (ephmax > epet_com) then
                  call dqg32(epet_com, ephmax, djpede, jpe(jz))
               else
                  jpe(jz) = 0.0_wp
               end if
            end if

            je(jz+1) = je0*sep*(1.0_wp + real(jz+1, wp)/tau)* &
                       (1.0_wp + sqrt(2.0_wp/(tau + real(2*jz+1, wp))))
            if (jz == 0) then
               ji(0) = ji0*(1.0_wp + sqrt(0.5_wp*PI/tau))
            else
               term = theta(jz)/tau
               if (term < ARGMAX) then
                  ji(jz) = ji0*(1.0_wp + 1.0_wp/sqrt(4.0_wp*tau + real(3*jz-1, wp)))**2 &
                           *exp(-term)
               else
                  ji(jz) = 0.0_wp
               end if
            end if
            fz(jz+1) = fz(jz)*(jpe(jz) + ji(jz))/je(jz+1)
            izmax = jz+1
            ssum  = ssum + fz(jz+1)
            if (fz(jz+1) < FMIN*ssum) exit
         end do
      end if

      ! ensure thresholds for state izmax
      call thresholds(icomp, izmax, a, ip(izmax), emin(izmax), epdt(izmax), epet(izmax))

      ! normalize fz, compute mode
      fac = 1.0_wp/ssum
      izmode = izmin
      do jz = -jzmx, izmin-1
         fz(jz) = 0.0_wp
      end do
      do jz = izmax+1, jzmx
         fz(jz) = 0.0_wp
      end do
      do jz = izmin, izmax
         fz(jz) = fac*fz(jz)
         if (fz(jz) > fz(izmode)) izmode = jz
      end do

      ! trim izmin/izmax around the mode
      if (izmin < izmode) then
         jzmin = izmin
         do jz = izmode, izmin, -1
            jzmin = jz
            if (fz(jz) < 1.0e-20_wp) exit
         end do
         izmin = jzmin
      end if
      if (izmax > izmode) then
         jzmax = izmax
         do jz = izmode, izmax
            jzmax = jz
            if (fz(jz) < 1.0e-20_wp) exit
         end do
         izmax = jzmax
      end if

      ! The HEAT / GPE collisional-cooling computation in the F77 source is
      ! not needed for the ionization fraction and is intentionally omitted.
      ! All quantities required by pah_ionfrac (fz, izmin, izmax) are set.
   end subroutine charge


   ! ======================================================================
   !  THRESHOLDS : photon-energy thresholds (Weingartner & Draine 1999)
   ! ======================================================================
   subroutine thresholds(icomp, iz, a, ip, emin, epdt, epet)
      integer,  intent(in)  :: icomp, iz
      real(wp), intent(in)  :: a
      real(wp), intent(out) :: ip, emin, epdt, epet
      real(wp) :: ebg, ipvb, w, x, z

      z = real(iz, wp)

      ! EMIN (composition-independent)
      if (iz >= -1) then
         emin = 0.0_wp
      else
         x    = (a/27.0e-8_wp)**0.75_wp
         emin = -(z + 1.0_wp)*(14.4e-8_wp/a)*x/(1.0_wp + x)
      end if

      if (icomp <= 2 .or. (icomp >= 11 .and. icomp <= 16)) then
         ! carbonaceous
         w   = 4.4_wp
         ebg = 0.0_wp
         ipvb = w + (14.4e-8_wp/a)*(z + 0.5_wp + (z + 2.0_wp)*0.3e-8_wp/a)
         if (iz < 0) then
            epet = max(ipvb + emin, 0.0_wp)
            ip   = w + (14.4e-8_wp/a)*(z + 0.5_wp - 4.0e-8_wp/(a + 7.0e-8_wp))
            epdt = max(ip + emin, 0.0_wp)
         else
            ip   = ipvb
            epet = ipvb
            epdt = epet
         end if
      else if (icomp == 3 .or. icomp == 7) then
         ! silicate
         w   = 8.0_wp
         ebg = 5.0_wp
         ipvb = w + (14.4e-8_wp/a)*(z + 0.5_wp + (z + 2.0_wp)*0.3e-8_wp/a)
         if (iz < 0) then
            ip   = w - ebg + (14.4e-8_wp/a)*(z + 0.5_wp)
            epdt = max(ip + emin, 0.0_wp)
            epet = max(ipvb + emin, 0.0_wp)
         else
            ip   = ipvb
            epet = ipvb
            epdt = epet
         end if
      else
         ip = 0.0_wp; epet = 0.0_wp; epdt = 0.0_wp
      end if
   end subroutine thresholds


   ! ======================================================================
   !  PEYIELD : photoelectric yield (Weingartner & Draine 1999)
   ! ======================================================================
   subroutine peyield(ephoton, icomp, arad, tgr, iz, emin, epet, peyld, peke)
      real(wp), intent(in)  :: ephoton, arad, tgr, emin, epet
      integer,  intent(in)  :: icomp, iz
      real(wp), intent(out) :: peyld, peke
      real(wp) :: alpha, aphys, beta, la, le, phi, qq, term, wbulk, y0, y1, y2, z
      real(wp) :: enim, wave_um

      z       = real(iz, wp)
      wave_um = 1.23984_wp/ephoton

      if (icomp <= 2 .or. (icomp >= 11 .and. icomp <= 16)) then
         wbulk = 4.4_wp
         aphys = arad
      else if (icomp == 3 .or. icomp == 7) then
         wbulk = 8.0_wp
         aphys = arad
      else
         peyld = 0.0_wp; peke = 0.0_wp
         return
      end if

      if (ephoton <= epet) then
         peyld = 0.0_wp
         peke  = 0.0_wp
         return
      end if

      ! Coulomb barrier outside double layer
      phi = 14.4e-8_wp*(z + 1.0_wp)/aphys

      if (iz < 0) then
         term = (ephoton - epet)/wbulk
      else
         term = (ephoton - epet + phi)/wbulk
      end if

      if (icomp <= 2 .or. (icomp >= 11 .and. icomp <= 16)) then
         term = term**5
         y0   = 9.0e-3_wp*term/(1.0_wp + 3.7e-2_wp*term)
      else
         y0   = 0.5_wp*term/(1.0_wp + 5.0_wp*term)
      end if

      ! Y2 and mean photoelectron KE
      if (iz < 0) then
         y2   = 1.0_wp
         peke = 0.5_wp*(ephoton - epet) + emin
      else
         qq   = phi/(ephoton - epet + phi)
         y2   = 1.0_wp - qq**2*(3.0_wp - 2.0_wp*qq)
         peke = 0.5_wp*(ephoton - epet + phi)*(1.0_wp + qq**3*(3.0_wp*qq - 4.0_wp))/ &
                (1.0_wp + qq**2*(2.0_wp*qq - 3.0_wp)) - phi
      end if

      ! Y1 small-particle enhancement.  For ICOMP=-1,11-16 use the
      ! 1/3-2/3 graphite Im(n) average exactly as the F77 INDEX calls.
      if (icomp == -1 .or. (icomp >= 11 .and. icomp <= 16)) then
         enim = graphite_imn(arad*1.0e4_wp, wave_um)   ! (Im n1 + 2 Im n2)/3
      else
         enim = graphite_imn(arad*1.0e4_wp, wave_um)    ! (silicate not exercised)
      end if

      la   = 9.8663e-6_wp/(ephoton*enim)
      le   = 1.0e-7_wp
      beta = aphys/la
      alpha = beta + aphys/le
      y1 = (beta/alpha)**2* &
           (exp(-alpha) - 1.0_wp + alpha*(1.0_wp - 0.5_wp*alpha))/ &
           (exp(-beta)  - 1.0_wp + beta *(1.0_wp - 0.5_wp*beta))

      peyld = y2*min(1.0_wp, y0*y1)
   end subroutine peyield


   ! ======================================================================
   !  DJPEDE / DGPEDE : the DQG32 integrands (read module COMMON state)
   ! ======================================================================
   function djpede(e) result(f)
      real(wp), intent(in) :: e
      real(wp) :: f
      integer  :: icomp, jz
      real(wp) :: eev, wavnum, en, qabs, peyld, peke, sigma, x
      real(wp) :: ar, tg, emin, epdt, epet, wave_um

      eev    = e
      wavnum = 8065.6_wp*e
      icomp  = icomp_com
      jz     = jz_com
      emin   = emin_com
      epdt   = epdt_com
      epet   = epet_com

      ! charged PAH/graphitic -> ICOMP=12 (only affects QCOMP/PEYIELD
      ! material; cross sections here are graphite-based, see deviations)
      if (icomp_com == 11 .and. jz /= 0) icomp = 12
      if (icomp_com == 13 .and. jz /= 0) icomp = 14
      if (icomp_com == 15 .and. jz /= 0) icomp = 16

      ar = ar_com
      tg = tg_com

      call radfld_mmp(rftyp_com, rfchi_com, rfr_v_com, rftau_com, &
                      rfion_com, rfpar_com, wavnum, en)

      wave_um = 1.0e4_wp/wavnum
      call qabs_pahgra(ar*1.0e4_wp, wave_um, icomp, qabs)

      call peyield(eev, icomp, ar, tg, jz, emin, epet, peyld, peke)

      sigma = 3.14159_wp*ar**2*qabs*peyld

      ! photodetachment of excess negative charge
      if (jz < 0) then
         if (e > epdt) then
            x = (e - epdt)/3.0_wp
            sigma = sigma - real(jz, wp)*1.219e-17_wp*x/(1.0_wp + x**2/3.0_wp)**2
         else
            sigma = 0.0_wp
         end if
      end if

      f = 3.946e23_wp*eev**2*en*sigma
   end function djpede


   function dgpede(e) result(f)
      real(wp), intent(in) :: e
      real(wp) :: f
      integer  :: icomp, jz
      real(wp) :: eev, wavnum, en, qabs, peyld, peke, sigke, x, eke
      real(wp) :: ar, tg, emin, epdt, epet, wave_um

      eev    = e
      wavnum = 8065.6_wp*e
      icomp  = icomp_com
      jz     = jz_com
      emin   = emin_com
      epdt   = epdt_com
      epet   = epet_com

      if (icomp == 11 .and. jz /= 0) icomp = 12
      if (icomp == 13 .and. jz /= 0) icomp = 14
      if (icomp == 15 .and. jz /= 0) icomp = 16

      ar = ar_com
      tg = tg_com

      call radfld_mmp(rftyp_com, rfchi_com, rfr_v_com, rftau_com, &
                      rfion_com, rfpar_com, wavnum, en)

      wave_um = 1.0e4_wp/wavnum
      call qabs_pahgra(ar*1.0e4_wp, wave_um, icomp, qabs)

      call peyield(eev, icomp, ar, tg, jz, emin, epet, peyld, peke)

      sigke = 3.14159_wp*ar**2*qabs*peyld*peke

      if (jz < 0) then
         if (e > epdt) then
            x   = (e - epdt)/3.0_wp
            eke = e - epdt + emin
            sigke = sigke - real(jz, wp)*1.219e-17_wp*x*eke/(1.0_wp + x**2/3.0_wp)**2
         else
            sigke = 0.0_wp
         end if
      end if

      f = 6.3221e11_wp*eev**2*en*sigke
   end function dgpede


   ! ======================================================================
   !  RADFLD (MODE=4) : Mathis, Mezger & Panagia radiation field
   !  Returns EN = photon occupation number (angle/pol averaged).
   ! ======================================================================
   subroutine radfld_mmp(mode, chi, r_v, tauv, ilycon, alpha, wavnum, en)
      integer,  intent(in)  :: mode, ilycon
      real(wp), intent(in)  :: chi, r_v, tauv, alpha, wavnum
      real(wp), intent(out) :: en
      real(wp) :: eev, fac, dust

      ! No dust attenuation (TAUV=0 for DISM) and only MODE=4 is exercised.
      dust = 1.0_wp
      en   = 0.0_wp

      if (mode == 4) then
         eev = 1.23984e-4_wp*wavnum
         fac = 0.8115_wp*chi*dust
         if (eev > 13.60_wp) then
            en = 0.0_wp
         else if (eev > 11.20_wp .and. eev <= 13.60_wp) then
            en = fac*4.7353e-2_wp*3.328e-9_wp*eev**(-8.4172_wp)
         else if (eev > 9.26_wp  .and. eev <= 11.20_wp) then
            en = fac*4.7353e-2_wp*8.463e-13_wp*eev**(-5.0_wp)
         else if (eev > 5.04_wp  .and. eev <= 9.26_wp) then
            en = fac*4.7353e-2_wp*2.055e-14_wp*eev**(-3.3322_wp)
         else if (eev <= 5.04_wp) then
            en = fac*(1.00e-14_wp/(exp(1.43877_wp*wavnum/7500.0_wp) - 1.0_wp) + &
                      1.65e-13_wp/(exp(1.43877_wp*wavnum/4000.0_wp) - 1.0_wp) + &
                      4.00e-13_wp/(exp(1.43877_wp*wavnum/3000.0_wp) - 1.0_wp))
         end if
      end if

      ! silence unused-dummy warnings (r_v, tauv, ilycon, alpha unused for MODE=4)
      if (r_v < 0.0_wp .or. tauv < 0.0_wp .or. ilycon < 0 .or. alpha < -1.0_wp) en = en
   end subroutine radfld_mmp


   ! ======================================================================
   !  DQG32 : 32-point Gaussian quadrature (IBM SSP); integrand via interface
   ! ======================================================================
   subroutine dqg32(xl, xu, fct, y)
      real(wp), intent(in)  :: xl, xu
      procedure(integrand)  :: fct
      real(wp), intent(out) :: y
      real(wp) :: a, b, c

      a = 0.5_wp*(xu + xl)
      b = xu - xl
      c = 0.49863193092474078_wp*b
      y = 0.35093050047350483e-2_wp*(fct(a+c) + fct(a-c))
      c = 0.49280575577263417_wp*b
      y = y + 0.8137197365452835e-2_wp*(fct(a+c) + fct(a-c))
      c = 0.48238112779375322_wp*b
      y = y + 0.12696032654631030e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.46745303796886984_wp*b
      y = y + 0.17136931456510717e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.44816057788302606_wp*b
      y = y + 0.21417949011113340e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.42468380686628499_wp*b
      y = y + 0.25499029631188088e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.39724189798397120_wp*b
      y = y + 0.29342046739267774e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.36609105937014484_wp*b
      y = y + 0.32911111388180923e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.33152213346510760_wp*b
      y = y + 0.36172897054424253e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.29385787862038116_wp*b
      y = y + 0.39096947893535153e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.25344995446611470_wp*b
      y = y + 0.41655962113473378e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.21067563806531767_wp*b
      y = y + 0.43826046502201906e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.16593430114106382_wp*b
      y = y + 0.45586939347881942e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.11964368112606854_wp*b
      y = y + 0.46922199540402283e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.7223598079139825e-1_wp*b
      y = y + 0.47819360039637430e-1_wp*(fct(a+c) + fct(a-c))
      c = 0.24153832843869158e-1_wp*b
      y = b*(y + 0.48270044257363900e-1_wp*(fct(a+c) + fct(a-c)))
   end subroutine dqg32


   ! ======================================================================
   !  Graphite optical-constant helpers (self-contained; mirror q_graphite)
   ! ======================================================================
   subroutine init_module()
      call load_tables()
      initialized = .true.
   end subroutine init_module


   subroutine load_tables()
      integer  :: i, u
      real(wp) :: ener, rn1, rk, e1, e2

      open(newunit=u, file=F_CPA, status='old', action='read')
      read(u, '(/)')
      do i = 1, NCPA
         read(u, *) ener, rn1, rk, e1, e2
         cpa_eV(i)   = ener
         cpa_eps1(i) = e1
         cpa_eps2(i) = e2
         cpa_wavl(i) = 1.23984_wp / ener   ! [um]
      end do
      close(u)

      open(newunit=u, file=F_CPE, status='old', action='read')
      read(u, '(/)')
      do i = 1, NCPE
         read(u, *) ener, rn1, rk, e1, e2
         cpe_eV(i)   = ener
         cpe_eps1(i) = e1
         cpe_eps2(i) = e2
         cpe_wavl(i) = 1.23984_wp / ener
      end do
      close(u)
   end subroutine load_tables


   subroutine free_diel(ener_eV, agrain, dtype, eps1, eps2)
      ! Free-electron (Drude) contribution to graphite dielectric function.
      real(wp),         intent(in)  :: ener_eV, agrain   ! agrain in cm
      character(len=*), intent(in)  :: dtype
      real(wp),         intent(out) :: eps1, eps2
      real(wp) :: e_plasma, tau_bulk, veff, tau, x, g

      if (dtype == 'CpeD03') then
         e_plasma = 0.285_wp*sqrt(1.0_wp - 6.24e-3_wp*T_GRAIN + 3.66e-5_wp*T_GRAIN*T_GRAIN)
         tau_bulk = 4.2e-11_wp/(1.0_wp + 0.322_wp*T_GRAIN + 1.30e-3_wp*T_GRAIN*T_GRAIN)
         veff     = 4.5e11_wp*sqrt(1.0_wp + T_GRAIN/255.0_wp)
      else
         e_plasma = 0.101_wp
         tau_bulk = 3.0e-14_wp
         veff     = 3.7e10_wp*sqrt(1.0_wp + T_GRAIN/255.0_wp)
      end if

      tau  = 1.0_wp/(1.0_wp/tau_bulk + veff/agrain)
      x    = ener_eV/e_plasma
      g    = 1.0_wp/(1.518e15_wp*e_plasma*tau)
      eps1 = -1.0_wp/(x*x + g*g)
      eps2 =  g/(x*(x*x + g*g))
   end subroutine free_diel


   subroutine eps_to_nk(eps1, eps2, n, k)
      ! eps1 is Re(eps)-1 (Draine convention).
      real(wp), intent(in)  :: eps1, eps2
      real(wp), intent(out) :: n, k
      real(wp) :: rr

      if ((eps1*eps1 + eps2*eps2) < 1.0e-6_wp) then
         n = 0.5_wp*eps1 - 0.125_wp*(eps1*eps1 - eps2*eps2) + 1.0_wp
         k = 0.5_wp*eps2 - 0.25_wp*eps1*eps2
      else if (eps2 < 1.0e-3_wp*abs(1.0_wp + eps1)) then
         n = sqrt(1.0_wp + eps1)*(1.0_wp + 0.125_wp*(eps2/(1.0_wp + eps1))**2)
         k = sqrt(1.0_wp + eps1)*0.5_wp*(eps2/(1.0_wp + eps1))
      else
         rr = sqrt((1.0_wp + eps1)**2 + eps2*eps2)
         n  = sqrt(0.5_wp*(rr + 1.0_wp + eps1))
         k  = sqrt(0.5_wp*(rr - 1.0_wp - eps1))
      end if
   end subroutine eps_to_nk


   subroutine build_nk(agrain)
      ! Refresh cached (n,k) tables for graphite-equivalent radius agrain [cm].
      real(wp), intent(in) :: agrain
      integer  :: i
      real(wp) :: e1f, e2f, e1, e2

      do i = 1, NCPA
         call free_diel(cpa_eV(i), agrain, 'CpaD03', e1f, e2f)
         e1 = cpa_eps1(i) + e1f
         e2 = cpa_eps2(i) + e2f
         call eps_to_nk(e1, e2, cpa_n(i), cpa_k(i))
      end do
      do i = 1, NCPE
         call free_diel(cpe_eV(i), agrain, 'CpeD03', e1f, e2f)
         e1 = cpe_eps1(i) + e1f
         e2 = cpe_eps2(i) + e2f
         call eps_to_nk(e1, e2, cpe_n(i), cpe_k(i))
      end do
      cached_a = agrain
   end subroutine build_nk


   function graphite_imn(a_um, lambda_um) result(enim)
      ! 1/3-2/3 averaged graphite Im(n): (Im n_parallel + 2 Im n_perp)/3.
      real(wp), intent(in) :: a_um, lambda_um
      real(wp) :: enim
      real(wp) :: a_cm, k_pa, k_pe

      if (.not. initialized) call init_module()
      a_cm = a_um*1.0e-4_wp
      if (a_cm /= cached_a) call build_nk(a_cm)

      call interp_lin(cpa_wavl, cpa_k, NCPA, lambda_um, k_pa)
      call interp_lin(cpe_wavl, cpe_k, NCPE, lambda_um, k_pe)
      enim = (k_pa + 2.0_wp*k_pe)/3.0_wp
   end function graphite_imn


   subroutine graphite_qabs(a_um, lambda_um, qabs)
      ! Random-orientation-averaged graphite Q_abs = (1/3)Q_par + (2/3)Q_perp.
      real(wp), intent(in)  :: a_um, lambda_um
      real(wp), intent(out) :: qabs
      real(wp) :: a_cm, x, n_pa, k_pa, n_pe, k_pe
      real(wp) :: qext, qsca, qabs_pa, qabs_pe, albe, gsca

      if (.not. initialized) call init_module()
      a_cm = a_um*1.0e-4_wp
      if (a_cm /= cached_a) call build_nk(a_cm)

      call interp_lin(cpa_wavl, cpa_n, NCPA, lambda_um, n_pa)
      call interp_lin(cpa_wavl, cpa_k, NCPA, lambda_um, k_pa)
      call interp_lin(cpe_wavl, cpe_n, NCPE, lambda_um, n_pe)
      call interp_lin(cpe_wavl, cpe_k, NCPE, lambda_um, k_pe)

      x = 2.0_wp*PI*a_um/lambda_um
      call mie_q(n_pa, k_pa, x, qext, qsca, qabs_pa, albe, gsca)
      call mie_q(n_pe, k_pe, x, qext, qsca, qabs_pe, albe, gsca)
      qabs = qabs_pa/3.0_wp + 2.0_wp*qabs_pe/3.0_wp
   end subroutine graphite_qabs


   subroutine qabs_pahgra(a_um, lambda_um, icomp, qabs)
      ! QCOMP ICOMP=11-16 branch: PAH/graphitic Q_abs.
      !   lambda <= 0.058 um : pure 1/3-2/3 graphite (UV; PAH==graphite)
      !   lambda >  0.058 um : FGMIN*graphite + (1-FGMIN)*[FAC*PAH + (1-FAC)*gra]
      ! with FAC = (A_T/a)^3 for a > A_T, FAC = 1 for a <= A_T.
      real(wp), intent(in)  :: a_um, lambda_um
      integer,  intent(in)  :: icomp
      real(wp), intent(out) :: qabs
      real(wp) :: q_gra, q_pah, cabs_pah, xfrac, fac, a_cm
      real(wp) :: n_c, n_ring
      logical  :: ionized

      call graphite_qabs(a_um, lambda_um, q_gra)

      if (lambda_um <= 0.05800_wp) then
         qabs = q_gra
         return
      end if

      ! number of C atoms (graphite-equivalent radius), capped at a=0.15 um
      if (a_um < 0.15_wp) then
         n_c = nint(4.68e11_wp*a_um**3 + 0.5_wp)
      else
         n_c = nint(4.68e11_wp*0.15_wp**3 + 0.5_wp)
      end if
      n_c = max(n_c, 1.0_wp)

      if (n_c <= 40.0_wp) then
         n_ring = 0.3_wp*n_c
      else
         n_ring = 0.4_wp*n_c
      end if

      ! charged PAH/graphitic states (12,14,16) use the ionized edge
      ionized = (icomp == 12 .or. icomp == 14 .or. icomp == 16)

      call pah_cabs_uv(lambda_um, n_ring, ionized, cabs_pah)   ! [cm^2/C-atom]

      ! Q_pah = N_C * C_abs / (pi a^2)   [a in cm]
      a_cm  = a_um*1.0e-4_wp
      q_pah = n_c*cabs_pah/(PI*a_cm**2)

      if (a_um <= A_T) then
         xfrac = FGMIN
      else
         fac   = (A_T/a_um)**3
         xfrac = FGMIN + (1.0_wp - FGMIN)*(1.0_wp - fac)
      end if

      qabs = xfrac*q_gra + (1.0_wp - xfrac)*q_pah
   end subroutine qabs_pahgra


   subroutine pah_cabs_uv(lambda_um, n_ring, ionized, cabs)
      ! UV/visible part of the W&D(2001)/PAH_CRS_SCT cross section
      ! [cm^2 per C atom].  Only the UV continuum is needed for charging;
      ! the IR vibrational features (lambda > ~1 um) play no role since
      ! charging photons have E > 5 eV.
      real(wp), intent(in)  :: lambda_um, n_ring
      logical,  intent(in)  :: ionized
      real(wp), intent(out) :: cabs
      real(wp), parameter :: width1 = 2.70_wp, xpeak1 = 13.85_wp
      real(wp), parameter :: width2 = 1.0_wp,  xpeak2 = 4.6_wp
      real(wp), parameter :: c1 = 1.8687_wp, c2 = 0.1905_wp
      real(wp), parameter :: c3 = 7.850_wp,  c4 = 0.7743_wp
      real(wp) :: x, crs, xc, x2xc, cx2xc, cs33, cs_visual

      x = 1.0_wp/lambda_um   ! 1/micron

      if (x > 15.0_wp .and. x <= 17.25_wp) then
         crs = 126.0_wp - 6.4943_wp*x
      else if (x > 10.0_wp .and. x <= 15.0_wp) then
         crs = -3.0_wp + 1.35_wp*x + 18.8_wp*width1**2/ &
               ((x - xpeak1**2/x)**2 + width1**2)
      else if (x > 7.7_wp .and. x <= 10.0_wp) then
         crs = 66.302_wp - 24.367_wp*x + 2.9501_wp*x*x - 0.10569_wp*x**3
      else if (x > 5.9_wp .and. x <= 7.7_wp) then
         crs = c1 + c2*x + c3/((x - xpeak2**2/x)**2 + width2**2) + &
               c4*(0.5392_wp*(x - 5.9_wp)**2 + 0.05644_wp*(x - 5.9_wp)**3)
      else if (x >= 3.3_wp .and. x <= 5.9_wp) then
         crs = c1 + c2*x + c3/((x - xpeak2**2/x)**2 + width2**2)
      else if (x >= 0.0_wp .and. x < 3.3_wp) then
         if (ionized) then
            xc = 1.0e-4_wp*(22824.0_wp*n_ring**(-0.5_wp) + 8892.0_wp)
         else
            xc = 1.0e-4_wp*(38040.0_wp*n_ring**(-0.5_wp) + 10520.0_wp)
         end if
         x2xc  = x/xc
         cx2xc = atan(1000.0_wp*((x2xc - 1.0_wp)**3)/x2xc)/3.1416_wp + 0.5_wp
         cs33  = c1 + c2*3.3_wp + c3*width2**2/ &
                 ((3.3_wp - xpeak2**2/3.3_wp)**2 + width2**2)
         cs_visual = 4.563e-17_wp*10.0_wp**(-3.431_wp/3.3_wp)/1.0e-18_wp
         crs = 4.563e-17_wp*10.0_wp**(-3.431_wp/x)/1.0e-18_wp*cx2xc* &
               cs33/cs_visual
      else
         crs = 0.0_wp
      end if

      cabs = crs*1.0e-18_wp   ! cm^2 / C atom
   end subroutine pah_cabs_uv


   subroutine interp_lin(x, y, n, xnew, ynew)
      ! Linear interpolation; x may be ascending or descending (locate-based).
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: x(n), y(n), xnew
      real(wp), intent(out) :: ynew
      integer  :: jl, ju, jm
      logical  :: ascend

      ascend = (x(n) >= x(1))
      jl = 0
      ju = n + 1
      do while (ju - jl > 1)
         jm = (ju + jl)/2
         if (ascend .eqv. (xnew >= x(jm))) then
            jl = jm
         else
            ju = jm
         end if
      end do

      if (jl <= 0) then
         ynew = y(1)
      else if (jl >= n) then
         ynew = y(n)
      else
         if (ascend) then
            ynew = y(jl)   + (y(jl+1)-y(jl))*(xnew-x(jl))/(x(jl+1)-x(jl))
         else
            ynew = y(jl+1) + (y(jl)-y(jl+1))*(xnew-x(jl+1))/(x(jl)-x(jl+1))
         end if
      end if
   end subroutine interp_lin


   subroutine mie_q(refn, refk, x, qext, qsca, qabs, albe, gsca)
      ! Bohren-Huffman Mie (BHMIE), Q-only.  Mirrors astrodust/sed/src/mie.f90.
      real(wp), intent(in)  :: refn, refk, x
      real(wp), intent(out) :: qext, qsca, qabs, albe, gsca
      integer  :: n, nstop, nmx, nn
      real(wp) :: chi, chi0, chi1, en, fn, p, psi, psi0, psi1, xstop, ymod
      real(wp) :: amu, pii, pi0, pi1, tau
      complex(wp) :: dcxs1, an, an1, bn, bn1, refrl, xi, xi1, y
      complex(wp), allocatable :: d(:)

      if (x == 0.0_wp) then
         qext = 0.0_wp; qsca = 0.0_wp; qabs = 0.0_wp; albe = 0.0_wp; gsca = 0.0_wp
         return
      end if

      refrl = cmplx(refn, refk, kind=wp)
      y     = x*refrl
      ymod  = abs(y)

      xstop = x + 4.0_wp*x**0.3333_wp + 2.0_wp
      nmx   = nint(max(xstop, ymod)) + 15
      nstop = nint(xstop)
      allocate(d(nmx))

      amu = 1.0_wp; pi0 = 0.0_wp; pi1 = 1.0_wp
      dcxs1 = (0.0_wp, 0.0_wp)
      d(nmx) = (0.0_wp, 0.0_wp)
      nn = nmx - 1
      do n = 1, nn
         en = nmx - n + 1
         d(nmx-n) = (en/y) - (1.0_wp/(d(nmx-n+1) + en/y))
      end do

      psi0 = cos(x); psi1 = sin(x)
      chi0 = -sin(x); chi1 = cos(x)
      xi1  = cmplx(psi1, -chi1, kind=wp)
      qsca = 0.0_wp; gsca = 0.0_wp; p = -1.0_wp
      an = (0.0_wp,0.0_wp); bn = (0.0_wp,0.0_wp)
      an1 = (0.0_wp,0.0_wp); bn1 = (0.0_wp,0.0_wp)
      do n = 1, nstop
         en = n
         fn = (2.0_wp*en + 1.0_wp)/(en*(en + 1.0_wp))
         psi = (2.0_wp*en - 1.0_wp)*psi1/x - psi0
         chi = (2.0_wp*en - 1.0_wp)*chi1/x - chi0
         xi  = cmplx(psi, -chi, kind=wp)
         if (n > 1) then
            an1 = an; bn1 = bn
         end if
         an = (d(n)/refrl + en/x)*psi - psi1
         an = an/((d(n)/refrl + en/x)*xi - xi1)
         bn = (refrl*d(n) + en/x)*psi - psi1
         bn = bn/((refrl*d(n) + en/x)*xi - xi1)
         qsca = qsca + (2.0_wp*en + 1.0_wp)*(abs(an)**2 + abs(bn)**2)
         gsca = gsca + ((2.0_wp*en + 1.0_wp)/(en*(en + 1.0_wp)))* &
                       (real(an, wp)*real(bn, wp) + aimag(an)*aimag(bn))
         if (n > 1) then
            gsca = gsca + ((en - 1.0_wp)*(en + 1.0_wp)/en)* &
                          (real(an1, wp)*real(an, wp) + aimag(an1)*aimag(an) + &
                           real(bn1, wp)*real(bn, wp) + aimag(bn1)*aimag(bn))
         end if
         pii = pi1
         tau = en*amu*pii - (en + 1.0_wp)*pi0
         dcxs1 = dcxs1 + fn*(an*pii + bn*tau)
         p = -p
         psi0 = psi1; psi1 = psi
         chi0 = chi1; chi1 = chi
         xi1  = cmplx(psi1, -chi1, kind=wp)
         pi1  = ((2.0_wp*en + 1.0_wp)*amu*pii - (en + 1.0_wp)*pi0)/en
         pi0  = pii
      end do

      gsca = 2.0_wp*gsca/qsca
      qsca = (2.0_wp/(x*x))*qsca
      qext = (4.0_wp/(x*x))*real(dcxs1, wp)
      qabs = qext - qsca
      albe = qsca/qext
      deallocate(d)
   end subroutine mie_q

end module pah_ioniz_mod
! ----------------------------------------------------------------------
! Deviations from the F77 source:
!
!  * The 2798-line QCOMP package is not ported wholesale, but its
!    ICOMP=11/12 branch IS reproduced: the graphite 1/3-2/3 Mie Q_abs
!    (from the Draine 2003 dielectric with the size-dependent free-electron
!    Drude term, as in q_graphite.f90) is blended with the W&D(2001) PAH
!    UV cross section (PAH_CRS_SCT, UV branch only) via the same
!    FGMIN/(A_T/a)^3 transition QCOMP uses.  Only the UV part of the PAH
!    cross section is ported -- charging photons all have E > 5 eV
!    (lambda < 0.25 um), so the IR vibrational features are never sampled.
!    The neutral/ionized (ICOMP 11 vs 12) distinction only changes the
!    visible absorption edge (x < 3.3 um^-1), which is likewise outside
!    the charging band; it is still threaded through for completeness.
!
!  * The Y1 small-particle enhancement Im(n) (the explicit INDEX wiring)
!    uses the identical graphite tables, matching the F77 PEYIELD exactly.
!
!  * RADFLD only implements MODE=4 (MMP), the single mode reached from
!    PAH_CHARGING_DISM; dust attenuation (TAUV=0) is omitted.
!
!  * The HEAT / GPE collisional-cooling block of CHARGE is omitted: it does
!    not affect fz, izmin, izmax, hence not the ionization fraction.
! ----------------------------------------------------------------------
