program test_rayleigh_zero_scattering
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use tmatrix_kinds, only: wp
   use asymptotic_optics, only: rayleigh_limit
   implicit none
   real(wp), parameter :: pi = acos(-1.0_wp)
   real(wp), parameter :: a = 0.01_wp, lam = 1.0_wp
   real(wp), parameter :: shapes(3) = [1.0_wp, 1.5_wp, 0.67_wp]
   real(wp) :: qext, qsca, walb, asymm, al1(3), al2(3), al3(3), al4(3), be1(3), be2(3)
   real(wp) :: qe_ori(3), qa_ori(3), qs_ori(3), qr_ori(3)
   real(wp) :: x, contrast, expected_qsca
   integer :: j, lmax

   do j = 1, size(shapes)
      call rayleigh_limit(a, lam, 1.0_wp, 0.0_wp, shapes(j), qext, qsca, walb, asymm, &
           al1, al2, al3, al4, be1, be2, lmax, qe_ori, qa_ori, qs_ori, qr_ori)
      if (qext /= 0.0_wp .or. qsca /= 0.0_wp .or. walb /= 0.0_wp) &
         error stop 'vacuum scalar cross sections'
      if (asymm /= 0.0_wp .or. lmax /= 0) error stop 'vacuum angular convention'
      if (any(al1 /= 0.0_wp) .or. any(al2 /= 0.0_wp) .or. any(al3 /= 0.0_wp) .or. &
          any(al4 /= 0.0_wp) .or. any(be1 /= 0.0_wp) .or. any(be2 /= 0.0_wp)) &
         error stop 'vacuum matrix convention'
      if (any(qe_ori /= 0.0_wp) .or. any(qa_ori /= 0.0_wp) .or. &
          any(qs_ori /= 0.0_wp) .or. any(qr_ori /= 0.0_wp)) &
         error stop 'vacuum orientation cross sections'

      call rayleigh_limit(a, lam, 1.5_wp, 0.0_wp, shapes(j), qext, qsca, walb, asymm, &
           al1, al2, al3, al4, be1, be2, lmax, qe_ori, qa_ori, qs_ori, qr_ori)
      if (.not. all(ieee_is_finite([qext,qsca,walb,asymm,al1,al2,al3,al4,be1,be2]))) &
         error stop 'nonfinite lossless output'
      if (qext <= 0.0_wp .or. qext /= qsca .or. walb /= 1.0_wp .or. lmax /= 2) &
         error stop 'lossless scattering limit'
      if (abs((qe_ori(1)+qe_ori(2)+qe_ori(3))/3.0_wp-qext) > 1e-14_wp*qext) &
         error stop 'orientation extinction average'
      if (abs((qs_ori(1)+qs_ori(2)+qs_ori(3))/3.0_wp-qsca) > 1e-14_wp*qsca) &
         error stop 'orientation scattering average'
      if (any(qa_ori /= 0.0_wp)) error stop 'lossless absorption'

      call rayleigh_limit(a, lam, 1.5_wp, 0.2_wp, shapes(j), qext, qsca, walb, asymm, &
           al1, al2, al3, al4, be1, be2, lmax, qe_ori, qa_ori, qs_ori, qr_ori)
      if (.not. all(ieee_is_finite([qext,qsca,walb,asymm,al1,al2,al3,al4,be1,be2]))) &
         error stop 'nonfinite absorbing output'
      if (qext <= qsca .or. qsca <= 0.0_wp .or. walb <= 0.0_wp .or. walb >= 1.0_wp) &
         error stop 'absorbing scattering limit'
      if (any(qa_ori <= 0.0_wp)) error stop 'absorbing orientation cross sections'
   end do

   x = 2.0_wp*pi*a/lam
   contrast = (1.5_wp**2-1.0_wp)/(1.5_wp**2+2.0_wp)
   expected_qsca = (8.0_wp/3.0_wp)*x**4*contrast**2
   call rayleigh_limit(a, lam, 1.5_wp, 0.0_wp, 1.0_wp, qext, qsca, walb, asymm, &
        al1, al2, al3, al4, be1, be2, lmax)
   if (abs(qsca/expected_qsca-1.0_wp) > 1e-13_wp) &
      error stop 'spherical Rayleigh cross section'
   if (abs(al1(1)-1.0_wp)>1e-14_wp .or. abs(al1(3)-0.5_wp)>1e-14_wp) &
      error stop 'spherical Rayleigh phase moments'
   al2 = -1.0_wp
   lmax = -1
   call rayleigh_limit(a, lam, 1.0_wp, 0.0_wp, 1.0_wp, qext, qsca, walb, asymm, &
        al2=al2, lmax=lmax)
   if (any(al2 /= 0.0_wp) .or. lmax /= 0) &
      error stop 'independent optional moments at zero scattering'
   call rayleigh_limit(a, lam, 1.5_wp, 0.0_wp, 1.0_wp, qext, qsca, walb, asymm, &
        al2=al2, lmax=lmax)
   if (abs(al2(3)-3.0_wp)>1e-14_wp .or. lmax /= 2) &
      error stop 'independent optional moments for sphere'
   write(*,'(a)') 'vacuum, lossless, absorbing, and spheroid Rayleigh limits passed'
end program test_rayleigh_zero_scattering
