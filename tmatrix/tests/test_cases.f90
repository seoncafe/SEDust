module tmatrix_test_cases
   !! 120 direct-core cases: three supported shapes, five size scales, four
   !! wavelengths, and weak/strong absorption.  Inputs stay below the legacy
   !! storage limits so this suite establishes numerical, not overflow, data.
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter, public :: n_cases = 120
   real(real64), parameter, public :: tolerance = 1.0e-3_real64
   integer, parameter, public :: ndgs = 2
   real(real64), parameter :: radii(5) = [0.010_real64, 0.050_real64, 0.100_real64, &
                                          0.250_real64, 0.500_real64]
   real(real64), parameter :: wavelengths(4) = [0.250_real64, 0.550_real64, &
                                                 1.000_real64, 5.000_real64]
   public :: get_case

contains

   subroutine get_case(i, a, lambda, nr, ki, aspect_ratio, shape)
      integer, intent(in) :: i
      real(real64), intent(out) :: a, lambda, nr, ki, aspect_ratio
      integer, intent(out) :: shape
      integer :: within_shape, ia, il, iabs, ishape

      ishape = (i-1) / 40 + 1
      within_shape = mod(i-1, 40)
      ia   = mod(within_shape, 5) + 1
      il   = mod(within_shape/5, 4) + 1
      iabs = within_shape/20 + 1
      a = radii(ia)
      lambda = wavelengths(il)
      if (mod(ia+il, 2) == 0) then
         nr = 1.40_real64
      else
         nr = 1.70_real64
      end if
      if (iabs == 1) then
         ki = 0.001_real64
      else
         ki = 0.100_real64
      end if
      select case (ishape)
      case (1)
         shape = -1
         aspect_ratio = 1.4_real64
      case (2)
         shape = -2
         aspect_ratio = 1.4_real64
      case default
         shape = 2
         aspect_ratio = 0.10_real64
      end select
   end subroutine get_case

end module tmatrix_test_cases
