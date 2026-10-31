program test_cmb_field_balance
   ! Check the full-field convention used by the SED callers of calc_P.
   ! The finite-bin GD solver has energy closure, not exact Planck detailed
   ! balance. No background term belongs in cooling when it is in Jfield.
   use constants, only: wp
   use radfield, only: bbody
   use p_sub, only: calc_P, p_sub_setup, radiation_absorption
   implicit none
   integer, parameter :: nl = 240, nt = 100
   real(wp), parameter :: pi = 3.141592653589793238462643383279502884197_wp
   real(wp), parameter :: tol = 2.0e-12_wp
   real(wp) :: lam(nl), cabs(nl), field(nl), temp(nt), h(nt), kappb(nt)
   real(wp) :: p(nt), lnp(nt), absorption, emission, rel
   integer :: i, j, icase

   do i = 1, nl
      lam(i) = 30.0_wp * (1.0e4_wp/30.0_wp)**(real(i-1, wp)/real(nl-1, wp))
      cabs(i) = 1.0e-21_wp * (1000.0_wp/lam(i))
   end do
   do j = 1, nt
      temp(j) = 0.5_wp * (20.0_wp/0.5_wp)**(real(j-1, wp)/real(nt-1, wp))
      h(j) = 1.0e-17_wp * temp(j)**3
      do i = 1, nl
         field(i) = bbody(temp(j), lam(i))
      end do
      kappb(j) = radiation_absorption(lam, cabs, field)
   end do
   call p_sub_setup(lam)

   field = 0.0_wp
   call calc_P(lam, cabs, field, temp, kappb, h, p, lnp)
   if (p(1) /= 1.0_wp .or. any(p(2:) /= 0.0_wp)) &
      error stop 'zero field produced phantom heating'
   write(*,'(a)') 'zero field: P(first) = 1, no phantom heating'

   do icase = 1, 2
      do i = 1, nl
         if (icase == 1) then
            field(i) = bbody(2.725_wp, lam(i))
         else
            field(i) = 0.0_wp
            if (lam(i) < 400.0_wp) field(i) = 0.3_wp*bbody(8.0_wp, lam(i))
         end if
      end do
      call calc_P(lam, cabs, field, temp, kappb, h, p, lnp)
      absorption = 4.0_wp*pi*radiation_absorption(lam, cabs, field)
      emission = 4.0_wp*pi*sum(p(2:)*kappb(2:))
      rel = emission/((1.0_wp-p(nt))*absorption) - 1.0_wp
      if (abs(rel) > tol) error stop 'full-field energy closure failed'
      if (abs(sum(p)-1.0_wp) > tol) error stop 'probability normalization failed'
      if (p(nt) > 1.0e-8_wp) error stop 'top enthalpy bin is populated'
      if (icase == 1) then
         write(*,'(a,es12.4)') 'pure CMB: relative energy closure = ', rel
      else
         write(*,'(a,es12.4)') 'field without CMB: relative energy closure = ', rel
      end if
   end do
end program test_cmb_field_balance
