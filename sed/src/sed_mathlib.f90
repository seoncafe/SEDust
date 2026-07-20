module sed_mathlib
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64
  implicit none
  private sp,dp,wp, PI, TWOPI
  integer,  parameter :: sp = real32
  integer,  parameter :: dp = real64
  integer,  parameter :: wp = dp
  real(wp), parameter :: PI = 3.141592653589793238462643383279502884197_wp
  real(wp), parameter :: TWOPI = PI * 2.0_wp
  interface interp
     module procedure interp1
  end interface interp
contains
  !--------------------
  function first_location(arr) result(loc)
  implicit none
  logical, intent(in) :: arr(:)
  integer :: loc, n, i
  n   = size(arr)
  loc = 0
  do i=1,n
     if (arr(i)) then
        loc = i
        exit
     endif
  enddo
  end function first_location
  !--------------------
  function last_location(arr) result(loc)
  implicit none
  logical, intent(in) :: arr(:)
  integer :: loc, n, i
  n = size(arr)
  loc = 0
  do i=n,1,-1
     if (arr(i)) then
        loc = i
        exit
     endif
  enddo
  end function last_location
  !--------------------
  function minloc1(array) result(loc1)
  implicit none
  real(wp), intent(in) :: array(:)
  integer :: loc(1), loc1
  loc  = minloc(array)
  loc1 = loc(1)
  end function
  !--------------------
  function maxloc1(array) result(loc1)
  implicit none
  real(wp), intent(in) :: array(:)
  integer :: loc(1), loc1
  loc  = maxloc(array)
  loc1 = loc(1)
  end function
  !--------------------
  subroutine interp_eq(x,y,xnew,ynew,ix)
  implicit none
  real(wp),     intent(in)  :: x(:),y(:),xnew
  real(wp),     intent(out) :: ynew
  integer, optional, intent(out) :: ix
  !---------------
  ! local variable
  integer :: n, i
  real(kind=wp) :: dx

  n  = size(x)
  dx = x(2) - x(1)
  i  = int((xnew-x(1))/dx + 1)
  if (i <= 0) then
     ynew = y(1)
  else if (i >= n) then
     ynew = y(n)
  else
     ynew = y(i) + (y(i+1)-y(i))*(xnew-x(i))/dx
  endif
  if (present(ix)) then
     ix = i
  endif

  return
  end subroutine interp_eq
  !--------------------
  subroutine interp1(x,y,xnew,ynew,ix)
  implicit none
  real(wp),          intent(in)  :: x(:),y(:),xnew
  real(wp),          intent(out) :: ynew
  integer, optional, intent(out) :: ix
  !---------------
  ! local variable
  integer :: n, i
  logical :: ascend

  n      = size(x)
  ascend = (x(n) >= x(1))
  call locate(x,xnew,i)
  if (i == 0) then
     ynew = y(1)
  else if (i == n) then
     ynew = y(n)
  else
     if (ascend) then
        ynew = y(i)   + (y(i+1)-y(i))*(xnew-x(i))/(x(i+1)-x(i))
     else
        ynew = y(i+1) + (y(i)-y(i+1))*(xnew-x(i+1))/(x(i)-x(i+1))
     endif
  endif
  if (present(ix)) then
     ix = i
  endif
  return
  end subroutine interp1
  !--------------------
  subroutine interp0(x,y,xnew,ynew,ix)
  implicit none
  real(wp),          intent(in)  :: x(:),y(:),xnew
  real(wp),          intent(out) :: ynew
  integer, optional, intent(out) :: ix
  !---------------
  ! local variable
  integer :: n, i
  logical :: ascend

  n      = size(x)
  ascend = (x(n) >= x(1))
  call locate(x,xnew,i)
  if (i == 0) then
     ynew = 0.0_wp
  else if (i == n) then
     ynew = 0.0_wp
  else
     if (ascend) then
        ynew = y(i)   + (y(i+1)-y(i))*(xnew-x(i))/(x(i+1)-x(i))
     else
        ynew = y(i+1) + (y(i)-y(i+1))*(xnew-x(i+1))/(x(i)-x(i+1))
     endif
  endif
  if (present(ix)) then
     ix = i
  endif
  return
  end subroutine interp0
  !--------------------
  subroutine locate(xx,x,j)
  implicit none
  real(wp), intent(in)  :: xx(:), x
  integer,  intent(out) :: j
  ! Taken from Numerical Recipes
  ! Given an array xx(1:n), and given a value x, returns a value j such that xx(j) <= x < xx(j+1).
  ! xx(1:n) must be monotonic, either increasing or decreasing.
  ! j = 0 or j =n is returned to indicate that x is out of range.

  integer :: n,jl,jm,ju
  logical :: ascend

   n      = size(xx)
   ascend = (xx(n) >= xx(1))
   jl     = 0   ! Initialize lower
   ju     = n+1 ! and upper limits.
   do
      if (ju-jl <= 1) exit  ! If we are done,
      jm = (ju+jl)/2        ! compute a midpoint,
      if (ascend .eqv. (x >= xx(jm))) then
         jl = jm            ! and replace either the lower limit
      else
         ju = jm            ! or the upper limit, as appropriate.
      endif
   enddo
   if (x == xx(1)) then ! Then set the output
     j = 1
   else if (x == xx(n)) then
     j = n-1
   else
     j = jl
   endif
   return
  end subroutine locate
  !--------------------
  function logadd(lna,lnb) result(lnc)
  implicit none
  real(wp), intent(in) :: lna, lnb
  real(wp) :: lnc
  if (lna >= lnb) then
     lnc = lna + log(1.0_wp + exp((lnb-lna)))
  else
     lnc = lnb + log(1.0_wp + exp((lna-lnb)))
  endif
  return
  end function logadd
  !--------------------
  function logadd_many(lna) result(lnc)
  implicit none
  real(wp), intent(in) :: lna(:)
  real(wp) :: lnc
  integer  :: n, i
  n   = size(lna)
  lnc = lna(1)
  if (n == 1) return
  do i=2,n
     lnc = logadd(lnc,lna(i))
  enddo
  return
  end function logadd_many
  !--------------------
end module sed_mathlib
