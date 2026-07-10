module mie_mod
   ! Bohren-Huffman Mie scattering for a homogeneous isotropic sphere.
   ! Lightly wrapped from astrodust/tmatrix/test/mie.f90 (which was
   ! originally B. T. Draine's BHMIE; K.-I. Seon's 2012 modifications
   ! kept the Q-only output and allocatable D() array).
   use constants, only: wp
   implicit none
   private
   public :: mie

contains

   subroutine mie(refn, refk, x, qext, qsca, qabs, albe, gsca)
      real(wp), intent(in)  :: refn, refk, x
      real(wp), intent(out) :: qext, qsca, qabs, albe, gsca

      integer :: n, nstop, nmx, nn
      real(wp) :: chi, chi0, chi1, en, fn, p, psi, psi0, psi1, xstop, ymod
      real(wp) :: amu, pi, pi0, pi1, tau
      complex(wp) :: dcxs1
      complex(wp) :: an, an1, bn, bn1, refrl, xi, xi1, y
      complex(wp), allocatable :: d(:)

      if (x == 0.0_wp) then
         qext = 0.0_wp
         qsca = 0.0_wp
         qabs = 0.0_wp
         albe = 0.0_wp
         gsca = 0.0_wp
         return
      end if

      refrl = cmplx(refn, refk, kind=wp)
      y     = x*refrl
      ymod  = abs(y)

      xstop = x + 4.0_wp*x**0.3333_wp + 2.0_wp
      nmx   = nint(max(xstop, ymod)) + 15
      nstop = nint(xstop)

      allocate(d(nmx))

      amu  = 1.0_wp
      pi0  = 0.0_wp
      pi1  = 1.0_wp
      dcxs1 = (0.0_wp, 0.0_wp)

      d(nmx) = (0.0_wp, 0.0_wp)
      nn = nmx - 1
      do n = 1, nn
         en = nmx - n + 1
         d(nmx-n) = (en/y) - (1.0_wp/(d(nmx-n+1) + en/y))
      end do

      psi0 = cos(x)
      psi1 = sin(x)
      chi0 = -sin(x)
      chi1 =  cos(x)
      xi1  = cmplx(psi1, -chi1, kind=wp)
      qsca = 0.0_wp
      gsca = 0.0_wp
      p    = -1.0_wp
      ! Initialize an, bn to avoid use-before-set on the n=1 iteration when
      ! the (n > 1) branches read an1/bn1 — they are gated by (n > 1), but
      ! some compilers still warn. Set explicit zeros.
      an = (0.0_wp, 0.0_wp)
      bn = (0.0_wp, 0.0_wp)
      an1 = (0.0_wp, 0.0_wp)
      bn1 = (0.0_wp, 0.0_wp)
      do n = 1, nstop
         en = n
         fn = (2.0_wp*en + 1.0_wp) / (en*(en + 1.0_wp))
         psi = (2.0_wp*en - 1.0_wp)*psi1/x - psi0
         chi = (2.0_wp*en - 1.0_wp)*chi1/x - chi0
         xi  = cmplx(psi, -chi, kind=wp)

         if (n > 1) then
            an1 = an
            bn1 = bn
         end if
         an = (d(n)/refrl + en/x)*psi - psi1
         an = an / ((d(n)/refrl + en/x)*xi - xi1)
         bn = (refrl*d(n) + en/x)*psi - psi1
         bn = bn / ((refrl*d(n) + en/x)*xi - xi1)

         qsca = qsca + (2.0_wp*en + 1.0_wp) * (abs(an)**2 + abs(bn)**2)
         gsca = gsca + ((2.0_wp*en + 1.0_wp)/(en*(en + 1.0_wp))) * &
                       (real(an, kind=wp)*real(bn, kind=wp) + aimag(an)*aimag(bn))
         if (n > 1) then
            gsca = gsca + ((en - 1.0_wp)*(en + 1.0_wp)/en) * &
                          (real(an1, kind=wp)*real(an, kind=wp) + aimag(an1)*aimag(an) + &
                           real(bn1, kind=wp)*real(bn, kind=wp) + aimag(bn1)*aimag(bn))
         end if

         pi  = pi1
         tau = en*amu*pi - (en + 1.0_wp)*pi0
         dcxs1 = dcxs1 + fn*(an*pi + bn*tau)

         p    = -p
         psi0 = psi1
         psi1 = psi
         chi0 = chi1
         chi1 = chi
         xi1  = cmplx(psi1, -chi1, kind=wp)
         pi1  = ((2.0_wp*en + 1.0_wp)*amu*pi - (en + 1.0_wp)*pi0) / en
         pi0  = pi
      end do

      gsca = 2.0_wp*gsca / qsca
      qsca = (2.0_wp/(x*x)) * qsca
      qext = (4.0_wp/(x*x)) * real(dcxs1, kind=wp)
      qabs = qext - qsca
      albe = qsca / qext

      deallocate(d)
   end subroutine mie

end module mie_mod
