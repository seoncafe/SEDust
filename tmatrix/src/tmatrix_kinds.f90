module tmatrix_kinds
   !! Numeric kind owned by libtmatrix.  This deliberately does not import
   !! SEDust's constants module.
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   integer, parameter, public :: wp = real64
end module tmatrix_kinds
