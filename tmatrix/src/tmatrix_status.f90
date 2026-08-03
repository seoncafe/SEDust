module tmatrix_status
   implicit none
   private

   integer, parameter, public :: TMATRIX_SUCCESS                    = 0
   integer, parameter, public :: TMATRIX_ERR_INVALID_INPUT          = 100
   integer, parameter, public :: TMATRIX_ERR_NOT_INITIALIZED        = 101
   integer, parameter, public :: TMATRIX_ERR_ALLOCATION             = 102
   integer, parameter, public :: TMATRIX_ERR_LEGACY_CORE            = 200
   integer, parameter, public :: TMATRIX_ERR_LAPACK_FACTOR          = 208
   integer, parameter, public :: TMATRIX_ERR_LAPACK_INVERT          = 209
   integer, parameter, public :: TMATRIX_ERR_NUMERICAL_NONFINITE    = 300
end module tmatrix_status
