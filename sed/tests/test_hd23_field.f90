program test_hd23_field
   use constants, only: wp, pi, c, h_si, c_si, k_si
   use radfield, only: J_HD23, J_Mathis, bbody, use_mathis_corrected
   implicit none
   real(wp), parameter :: lam(14) = [0.08_wp, 0.0912_wp, 0.10_wp, 0.109999_wp, &
        0.110_wp, 0.12_wp, 0.133999_wp, 0.134_wp, 0.20_wp, 0.245999_wp, &
        0.246_wp, 0.25_wp, 1.0_wp, 1000.0_wp]
   real(wp), parameter :: u = 1.6_wp
   real(wp) :: field(size(lam)), zero_field(size(lam)), reference, stellar
   real(wp) :: mathis_lam(5), mathis_field(5), expected, w4000, tcmb
   integer :: i, convention

   call J_HD23(u, lam, field)
   call J_HD23(0.0_wp, lam, zero_field)
   do i = 1, size(lam)
      ! Piecewise, as in Draine's mmp_radiation_field.f: the MMP83 UV fits up
      ! to and including 0.246 um, the three diluted blackbodies above it.
      if (lam(i) <= 0.246_wp) then
         stellar = uv_intensity(lam(i))
      else
         stellar = 7.0e-13_wp*planck(3000.0_wp, lam(i)) &
                 + 1.65e-13_wp*planck(4000.0_wp, lam(i)) &
                 + 1.0e-14_wp*planck(7500.0_wp, lam(i))
      end if
      reference = u*stellar + planck(2.725_wp, lam(i))
      call check_close(field(i), reference, 'HD23 field')
      call check_close(zero_field(i), planck(2.725_wp, lam(i)), 'unscaled CMB')
      if (lam(i) <= 1.0_wp) &
         call check_close(field(i)-zero_field(i), u*stellar, 'U scaling')
   end do

   mathis_lam = [0.08_wp, 0.10_wp, 0.20_wp, 0.5_wp, 1000.0_wp]
   do convention = 1, 2
      use_mathis_corrected = convention == 1
      if (use_mathis_corrected) then
         w4000 = 1.65e-13_wp
         tcmb = 2.725_wp
      else
         w4000 = 1.0e-13_wp
         tcmb = 2.9_wp
      end if
      call J_Mathis(u, mathis_lam, mathis_field)
      do i = 1, size(mathis_lam)
         if (mathis_lam(i) < 0.0912_wp) then
            expected = 0.0_wp
         else if (mathis_lam(i) < 0.110_wp) then
            expected = 3069.0_wp*mathis_lam(i)**3.4172_wp
         else if (mathis_lam(i) < 0.134_wp) then
            expected = 1.627_wp
         else if (mathis_lam(i) < 0.250_wp) then
            expected = 0.0566_wp*mathis_lam(i)**(-1.6678_wp)
         else
            expected = 1.0e-14_wp*bbody(7500.0_wp, mathis_lam(i)) &
                     + w4000*bbody(4000.0_wp, mathis_lam(i)) &
                     + 4.0e-13_wp*bbody(3000.0_wp, mathis_lam(i))
         end if
         expected = u*expected + bbody(tcmb, mathis_lam(i))
         call check_close(mathis_field(i), expected, 'legacy Mathis')
      end do
   end do
   use_mathis_corrected = .true.
   write(*,'(a)') 'HD23 boundaries, CMB, U scaling, and both Mathis conventions passed'

contains

   function uv_intensity(wavelength) result(value)
      real(wp), intent(in) :: wavelength
      real(wp) :: value, lambda_u
      lambda_u = 0.0_wp
      if (wavelength >= 0.0912_wp .and. wavelength < 0.110_wp) then
         lambda_u = 1.287e-9_wp*wavelength**4.4172_wp
      else if (wavelength >= 0.110_wp .and. wavelength < 0.134_wp) then
         lambda_u = 6.825e-13_wp*wavelength
      else if (wavelength >= 0.134_wp .and. wavelength <= 0.246_wp) then
         lambda_u = 2.373e-14_wp*wavelength**(-0.6678_wp)
      end if
      value = 0.0_wp
      if (wavelength > 0.0_wp) &
         value = c/(4.0_wp*pi)*1.0e3_wp*lambda_u/wavelength
   end function uv_intensity

   function planck(temperature, wavelength) result(value)
      real(wp), intent(in) :: temperature, wavelength
      real(wp) :: value, length_m, exponent
      length_m = wavelength*1.0e-6_wp
      exponent = h_si*c_si/(length_m*k_si*temperature)
      value = 0.0_wp
      if (exponent < 700.0_wp) &
         value = (2.0_wp*h_si*c_si**2/length_m**5)/(exp(exponent)-1.0_wp)
   end function planck

   subroutine check_close(actual, target, label)
      real(wp), intent(in) :: actual, target
      character(*), intent(in) :: label
      real(wp) :: tolerance
      tolerance = 2.0e-11_wp*max(abs(target), 1.0e-100_wp)
      if (abs(actual-target) > tolerance) then
         write(*,'(a,2es24.15)') label//': actual and expected: ', actual, target
         error stop 'field comparison failed'
      end if
   end subroutine check_close
end program test_hd23_field
