! 2013-05-01, Kwang-Il Seon
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
module enthalpy
! module to calculate enthalpies for graphite, silicate, and PAH.
!      Debye function (or integral) of order n = 2,3
use constants, only: wp
implicit none
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Enthalpies for silicate, graphite grains, and PAH
! Draine & Li (2001, ApJ, 551, 807)
! input: T in Kelvin
! ouput: erg
   function enthalpy_DL01(T,radius,dtype) result(U)
   implicit none
   real(kind=wp) :: U
   real(kind=wp), intent(in) :: T
   real(kind=wp), intent(in) :: radius
   character(len=4), intent(in) :: dtype
! Local variables
   real(kind=wp), parameter :: Top = 863.0_wp, Tip = 2504.0_wp, T2 = 500.0_wp, T3 = 1500.0_wp
   real(kind=wp), parameter :: kB  = 1.3806488e-16_wp ! erg K^-1
! C2 = h x c / k, lambdaj(1:3) in cm
! Note C2, and hwc should be double precision to avoid NaN.
   real(kind=wp), parameter :: C2 = 1.4387687_wp
   real(kind=wp), parameter :: lambdaj(3) = (/11.3_wp, 8.6_wp, 3.3_wp/) * 1e-4_wp
   real(kind=wp), parameter :: hwc(3) = C2/lambdaj(:)

   real(kind=wp) :: natom
   real(kind=wp) :: H_C, UCH

   select case (trim(dtype))
   case ('Car0', 'Car1')
      natom = 4.6820810e11_wp*radius**3
      U     = (natom-2._wp)*kB*T*(debye2(real(Top/T,8)) + 2.0_wp*debye2(real(Tip/T,8)))
   case('Sil')
      natom = 7d0*(5.1040126e10_wp*radius**3)
      U     = (natom-2._wp)*kB*T*(2.0_wp*debye2(real(T2/T,8)) + debye3(real(T3/T,8)))
   case default
      print*, '! enthalpy_DL01: calculations only available for Gra, PAH, or Sil'
      stop 1
   end select

   ! DL01-original C-H mode size threshold (natom <= 5.75e4, a <~ 50 AA).
   ! (Previously modified to 1.00e4; reverted to the Draine & Li 2001 value.)
   if ((dtype == 'Car0' .or. dtype == 'Car1') .and. natom <= 5.75e4_wp) then
      ! H/C = the hydrogen to carbon ratio
      ! NH = number of H
      if (natom <= 25) then
         H_C = 0.5_wp
         !NH = int(0.5*natom+0.5)
      elseif (natom < 100) then
         H_C = 0.5_wp/sqrt(real(natom,wp)/25._wp)
         !NH = int(2.5*sqrt(real(natom))+0.5)
      else
         H_C = 0.25_wp
         !NH = int(0.25*natom+0.5)
      endif
      ! C-H mode
      UCH = (H_C*natom) * kB * sum(hwc(:)/(exp(hwc(:)/T)-1.0_wp))
      !UCH = NH * kB * sum(hwc(:)/(exp(hwc(:)/T)-1.0))
      U = U + UCH
   endif
   end function enthalpy_DL01
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   FUNCTION cheval(n, a, t) RESULT(fn_val)
!   This function evaluates a Chebyshev series, using the Clenshaw method
!   with Reinsch modification, as analysed in the paper by Oliver.
!
!   INPUT PARAMETERS
!       N - INTEGER - The no. of terms in the sequence
!       A - REAL (wp) ARRAY, dimension 0 to N - The coefficients of
!           the Chebyshev series
!       T - REAL (wp) - The value at which the series is to be evaluated
!
!   REFERENCES
!        "An error analysis of the modified Clenshaw method for
!         evaluating Chebyshev and Fourier series" J. Oliver,
!         J.I.M.A., vol. 20, 1977, pp379-391
!
! MACHINE-DEPENDENT CONSTANTS: NONE
! INTRINSIC FUNCTIONS USED;
!    ABS
! AUTHOR:  Dr. Allan J. MacLeod,
!          Dept. of Mathematics and Statistics,
!          University of Paisley ,
!          High St.,
!          PAISLEY,
!          SCOTLAND
!
! LATEST MODIFICATION:   21 December , 1992
   
   INTEGER, INTENT(IN)    :: n
   REAL(kind=wp), INTENT(IN)  :: a(0:n)
   REAL(kind=wp), INTENT(IN)  :: t
   REAL(kind=wp)              :: fn_val
   
   INTEGER    :: i
   REAL(kind=wp)  :: d1, d2, tt, u0, u1, u2
   REAL(kind=wp), PARAMETER  :: zero = 0.0_8, half = 0.5_8, test = 0.6_8,  &
                               two = 2.0_8
   
   u1 = zero
   ! Init d2/u2 so the (compiler-visible) n<0 path cannot read them; the
   ! actual series always has n>=0 and overwrites them in the loops below.
   u2 = zero
   d2 = zero

   !   If ABS ( T )  < 0.6 use the standard Clenshaw method
   IF (ABS(t) < test) THEN
     u0 = zero
     tt = t + t
     DO  i = n, 0, -1
       u2 = u1
       u1 = u0
       u0 = tt * u1 + a(i) - u2
     END DO
     fn_val = (u0-u2) / two
   ELSE
   !   If ABS ( T )  > =  0.6 use the Reinsch modification
     d1 = zero
   !   T > =  0.6 code
     IF (t > zero) THEN
       tt = (t-half) - half
       tt = tt + tt
       DO  i = n, 0, -1
         d2 = d1
         u2 = u1
         d1 = tt * u2 + a(i) + d2
         u1 = d1 + u2
       END DO
       fn_val = (d1+d2) / two
     ELSE
   !   T < =  -0.6 code
       tt = (t+half) + half
       tt = tt + tt
       DO  i = n, 0, -1
         d2 = d1
         u2 = u1
         d1 = tt * u2 + a(i) - d2
         u1 = d1 - u2
       END DO
       fn_val = (d1-d2) / two
     END IF
   END IF
   RETURN
   END FUNCTION cheval
   function debye2 (xvalue)
!*****************************************************************************80
!! DEBYE2 calculates the Debye function of order 2.
!
!  Discussion:
!    The function is defined by:
!       DEBYE2(x) = 2 / x^2 * Integral ( 0 <= t <= x ) t^2 / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!    24 August 2004
!
!  Author:
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!    Input,  real(kind = 8) XVALUE, the argument of the function.
!    Output, real(kind = 8) DEBYE2, the value of the function.
!
     implicit none
   
!  real(kind = 8) cheval
     real(kind = 8) debye2
     real(kind = 8), parameter :: eight = 8.0D+00
     real(kind = 8), parameter :: four = 4.0D+00
     real(kind = 8), parameter :: half = 0.5D+00
     integer i
     integer nexp
     integer, parameter :: nterms = 17
     real(kind = 8), parameter :: one = 1.0D+00
     real(kind = 8), parameter :: three = 3.0D+00
     real(kind = 8), parameter :: two = 2.0D+00
     real(kind = 8) x
     real(kind = 8) xvalue
     real(kind = 8), parameter :: zero = 0.0D+00
   
     real(kind = 8) adeb2(0:18),debinf,expmx, &
          rk,sum1,t,twent4,xk,xlim1, &
          xlim2,xlow,xupper
     data twent4/24.0d0/
     data debinf/4.80822761263837714160d0/
     data adeb2/2.59438102325707702826d0, &
                0.28633572045307198337d0, &
               -0.1020626561580467129d-1, &
                0.60491097753468435d-3, &
               -0.4052576589502104d-4, &
                0.286338263288107d-5, &
               -0.20863943030651d-6, &
                0.1552378758264d-7, &
               -0.117312800866d-8, &
                0.8973585888d-10, &
               -0.693176137d-11, &
                0.53980568d-12, &
               -0.4232405d-13, &
                0.333778d-14, &
               -0.26455d-15, &
                0.2106d-16, &
               -0.168d-17, &
                0.13d-18, &
               -0.1d-19/
   !
   !   Machine-dependent constants
   !
     data xlow,xupper/0.298023d-7,35.35051d0/
     data xlim1,xlim2/708.39642d0,2.1572317d154/
   !
     x = xvalue
   
     if (x < zero) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'DEBYE2 - Fatal error!'
       write ( *, '(a)' ) '  Argument X < 0.'
       debye2 = zero
     else if (x < xlow) then
       debye2 = ((x - eight) * x + twent4) / twent4
     else if (x <= four) then
       t      = ((x * x / eight) - half) - half
       debye2 = cheval(nterms, adeb2, t) - x / three
     else if (x <= xupper) then
       expmx = exp(-x)
       sum1  = zero
       rk    = aint(xlim1/x)
       nexp  = int(rk)
       xk    = rk * x
       do i = nexp, 1, -1
         t    =  (one + two/xk + two/(xk * xk)) / rk
         sum1 = sum1 * expmx + t
         rk   = rk - one
         xk   = xk - x
       end do
       debye2 = debinf / (x * x) - two * sum1 * expmx
     else if (x < xlim1) then
       expmx  = exp (-x)
       sum1   = ((x + two) * x + two) / (x * x)
       debye2 = debinf/(x * x) - two * sum1 * expmx
     else if (x <= xlim2) then
       debye2 = debinf / (x * x)
     else
       debye2 = zero
     end if
   
     return
   end function debye2
   function debye3 (xvalue)
!*****************************************************************************80
!! DEBYE3 calculates the Debye function of order 3.
!
!  Discussion:
!    The function is defined by:
!      DEBYE3(x) = 3 / x^3 * Integral ( 0 <= t <= x ) t^3 / ( exp ( t ) - 1 ) dt
!
!    The code uses Chebyshev series whose coefficients
!    are given to 20 decimal places.
!
!    This subroutine is set up to work on IEEE machines.
!
!  Modified:
!    07 August 2004
!
!  Author:
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
!  Reference:
!    Allan McLeod,
!    Algorithm 757, MISCFUN: A software package to compute uncommon
!      special functions,
!    ACM Transactions on Mathematical Software,
!    Volume 22, Number 3, September 1996, pages 288-301.
!
!  Parameters:
!    Input,  real(kind = 8) XVALUE, the argument of the function.
!    Output, real(kind = 8) DEBYE3, the value of the function.
!
     implicit none
   
!  real (kind = 8) cheval
     real (kind = 8) debye3
     real (kind = 8), parameter :: eight = 8.0D+00
     real (kind = 8), parameter :: four = 4.0D+00
     real (kind = 8), parameter :: half = 0.5D+00
     integer i
     integer nexp
     integer, parameter :: nterms = 16
     real (kind = 8), parameter :: one = 1.0D+00
     real (kind = 8), parameter :: six = 6.0D+00
     real (kind = 8), parameter :: three = 3.0D+00
     real (kind = 8) x
     real (kind = 8) xvalue
     real (kind = 8), parameter :: zero = 0.0D+00
   
     real (kind = 8) adeb3(0:18),debinf,expmx, &
          pt375,rk,sevp5,sum1,t,twenty, &
          xk,xki,xlim1,xlim2,xlow,xupper
     data pt375/0.375d0/
     data sevp5,twenty/7.5d0 , 20.0d0/
     data debinf/0.51329911273421675946d-1/
     data adeb3/2.70773706832744094526d0, &
                0.34006813521109175100d0, &
               -0.1294515018444086863d-1, &
                0.79637553801738164d-3, &
               -0.5463600095908238d-4, &
                0.392430195988049d-5, &
               -0.28940328235386d-6, &
                0.2173176139625d-7, &
               -0.165420999498d-8, &
                0.12727961892d-9, &
               -0.987963459d-11, &
                0.77250740d-12, &
               -0.6077972d-13, &
                0.480759d-14, &
               -0.38204d-15, &
                0.3048d-16, &
               -0.244d-17, &
                0.20d-18, &
               -0.2d-19/
   !
   !   Machine-dependent constants
   !
     data xlow,xupper/0.298023d-7,35.35051d0/
     data xlim1,xlim2/708.39642d0,0.9487163d103/
   !
     x = xvalue
   
     if (x < zero) then
       write (*, '(a)') ' '
       write (*, '(a)') 'DEBYE3 - Fatal error!'
       write (*, '(a)') '  Argument X < 0.'
       debye3 = zero
       return
     end if
   
     if (x < xlow) then
       debye3 = ( ( x - sevp5 ) * x + twenty ) / twenty
     else if (x <= 4) then
       t      = ((x * x / eight) - half) - half
       debye3 = cheval(nterms, adeb3, t) - pt375 * x
     else
   !
   !   Code for x > 4.0
   !
        if (xlim2 < x) then
           debye3 = zero
        else
           debye3 = one / (debinf * x * x * x)
           if (x < xlim1) then
              expmx = exp (-x)
              if (xupper < x) then
                 sum1 = (((x + three) * x + six) * x + six) / (x * x * x)
              else
                 sum1 = zero
                 rk   = aint(xlim1 / x)
                 nexp = int(rk)
                 xk   = rk * x
                 do i = nexp, 1, -1
                    xki  = one / xk
                    t    =  (((six * xki + six) * xki + three) * xki + one) / rk
                    sum1 = sum1 * expmx + t
                    rk   = rk - one
                    xk   = xk - x
                 end do
              end if
              debye3 = debye3 - three * sum1 * expmx
           end if
        end if
     end if
   
     return
   end function debye3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! calculates the enthalpy for a spherical silicate or graphite grain
!       of radius a, temperature T, containing natom atoms
! Guhathakurta & Draine (1989)
!
! Units:
!     enthalpy : ergs
!     T        : Kelvin
!     radius   : micron
   function enthalpy_GD89(T,natom,radius,dtype) result(U)
   implicit none
   real(kind=wp) :: U
   real(kind=wp), intent(in) :: T,radius
!   integer(kind=8), intent(in) :: natom
   real(kind=wp), intent(in) :: natom
   character(len=3), intent(in) :: dtype

! local variables
! volume: grain volume = 4pi/3 * radius^3 (cm^3)
   real(kind=wp) :: volume

   select case (dtype)
   ! Silicate (Leger, Jura, & Omont 1985; G&D1989)
   case ('Sil')
      volume = 4.18878_wp*1.e-12_wp * radius**3
      if (T <= 50._wp) then
         U = U1(T)
      else if (T <= 150._wp) then
         U = U1(50._wp) + U2(T)
      else if (T <= 500._wp) then
         U = U1(50._wp) + U2(150._wp) + U3(T)
      else
         U = U1(50._wp) + U2(150._wp) + U3(500._wp) + U4(T)
      end if
      U = (1.0_wp-(2.0_wp/real(natom,8))) * volume * U
   ! Graphite (Chase et al. 1985; G&D1989)
   case('Gra')
      U = (4.15e-22_wp*T**3.3_wp)/(1.0_wp+6.51e-3_wp*T+1.5e-6_wp*T*T+8.3e-7_wp*T**2.3_wp)
      U = (real(natom,8)-2.0_wp) * U
   case default
      print*, '! enthalpy_GD89: calculations only available for Gra & Sil identifiers'
      stop 1
   end select

   contains
   ! Piecewise silicate enthalpy segments (G&D1989 fits)
   pure function U1(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = (1.4e3_wp/3.0_wp)  * T**3
   end function U1
   pure function U2(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = (2.2e4_wp/2.3_wp)  * (T**2.3_wp  - 50._wp**2.3_wp)
   end function U2
   pure function U3(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = (4.8e5_wp/1.68_wp) * (T**1.68_wp - 150._wp**1.68_wp)
   end function U3
   pure function U4(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = 3.41e7_wp * (T - 500.0_wp)
   end function U4
   end function enthalpy_GD89
!------------------------------------
! Enthalpy for graphite given by Dwek et al. (1986, ApJ, 302, 363)
!
   function enthalpy_Dwek(T,natom,radius) result(U)
   implicit none
   real(kind=wp) :: U
   real(kind=wp), intent(in) :: T,radius
!   integer(kind=8), intent(in) :: natom
   real(kind=wp), intent(in) :: natom
   real(kind=wp) :: volume

   if (T <= 60._wp) then
      U = U1(T)
   else if (T <=100._wp) then
      U = U1(60._wp) + U2(T)
   else if (T <= 470._wp) then
      U = U1(60._wp) + U2(100._wp) + U3(T)
   else if (T <=1070._wp) then
      U = U1(60._wp) + U2(100._wp) + U3(470._wp) + U4(T)
   else
      U = U1(60._wp) + U2(100._wp) + U3(470._wp) + U4(1070._wp) + U5(T)
   end if
   volume = 4.18878_wp*1.e-12_wp * radius**3
   U      = (1.0_wp-(2.0_wp/real(natom,8))) * volume * U

   contains
   ! Piecewise graphite enthalpy segments (Dwek 1986 fits)
   pure function U1(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = (3.84e2_wp/3.00_wp) * T**3
   end function U1
   pure function U2(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = (2.32e3_wp/2.56_wp) * (T**2.56_wp -  60._wp**2.56_wp)
   end function U2
   pure function U3(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = (5.61e3_wp/2.37_wp) * (T**2.37_wp - 100._wp**2.37_wp)
   end function U3
   pure function U4(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = (7.74e5_wp/1.57_wp) * (T**1.57_wp - 470._wp**1.57_wp)
   end function U4
   pure function U5(T) result(res)
      real(kind=wp), intent(in) :: T
      real(kind=wp) :: res
      res = 4.14e7_wp * (T - 1070._wp)
   end function U5
   end function enthalpy_Dwek
end module enthalpy
