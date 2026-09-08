module effective_medium
   ! Effective dielectric function of a two-component composite grain.
   !
   ! maxwell_garnett is the Maxwell Garnett mixing rule for a matrix of
   ! dielectric function eps_m holding a dilute dispersion of spherical
   ! inclusions of eps_i occupying a volume fraction f.  It follows from
   ! applying the Clausius-Mossotti relation to the polarizability of a small
   ! sphere embedded in the matrix, so it inherits that derivation's domain of
   ! validity:
   !
   !   * the inclusions are small compared with the wavelength inside the
   !     matrix, so each responds as a static dipole;
   !   * they are far enough apart that each sees only the average field,
   !     which restricts f to small values -- the rule is asymmetric in the
   !     two components and is the recommended one when one of them is
   !     dilute, whereas Bruggeman is the symmetric rule used when the two
   !     volumes are comparable (Bohren & Huffman 1983, Sect. 8.5).
   !
   ! In the form Bohren & Huffman give it,
   !
   !   (eps_eff - eps_m) / (eps_eff + 2 eps_m)
   !                             = f (eps_i - eps_m) / (eps_i + 2 eps_m),
   !
   ! whose solution is the expression coded below,
   !
   !   eps_eff = eps_m [ 1 + 3 f (eps_i - eps_m)
   !                          / (eps_i + 2 eps_m - f (eps_i - eps_m)) ].
   !
   ! Guillet et al. (2018, A&A 610, A16) Eq. (28) is this same expression
   ! written for the refractive index, m_eff = sqrt(eps_eff) with eps = m^2,
   ! and is what their Model D uses for the astrosilicate matrix carrying 6%
   ! by volume of amorphous-carbon inclusions.
   !
   ! Limits the coded form reproduces exactly: f = 0 gives eps_m, f = 1 gives
   ! eps_i, and eps_i = eps_m gives eps_m for every f.  maxwell_garnett_m
   ! wraps the same rule in the refractive index, m = n + i k with k >= 0;
   ! the principal square root is the physical branch because a passive
   ! medium has Im(eps) >= 0, which places sqrt(eps) in the first quadrant.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   implicit none
   private
   public :: maxwell_garnett, maxwell_garnett_m

   integer, parameter :: wp = real64

contains

   pure function maxwell_garnett(eps_matrix, eps_inclusion, f_inclusion) &
        result(eps_eff)
      !! Effective dielectric function of the composite.
      complex(wp), intent(in) :: eps_matrix, eps_inclusion
      real(wp),    intent(in) :: f_inclusion       ! volume fraction, 0 <= f <= 1
      complex(wp) :: eps_eff
      complex(wp) :: d

      d = eps_inclusion - eps_matrix
      eps_eff = eps_matrix * (1.0_wp + 3.0_wp*f_inclusion*d &
                / (eps_inclusion + 2.0_wp*eps_matrix - f_inclusion*d))
   end function maxwell_garnett


   subroutine maxwell_garnett_m(nr_matrix, ki_matrix, nr_inclusion, ki_inclusion, &
                                f_inclusion, nr_eff, ki_eff)
      !! The same rule expressed in the refractive index, m = n + i k.
      real(wp), intent(in)  :: nr_matrix, ki_matrix
      real(wp), intent(in)  :: nr_inclusion, ki_inclusion
      real(wp), intent(in)  :: f_inclusion
      real(wp), intent(out) :: nr_eff, ki_eff
      complex(wp) :: eps_m, eps_i, eps_eff, m_eff

      if (f_inclusion < 0.0_wp .or. f_inclusion > 1.0_wp) then
         write(error_unit,'(a,es12.4)') &
            'maxwell_garnett_m: volume fraction outside [0,1]: ', f_inclusion
         stop 1
      end if
      eps_m = cmplx(nr_matrix,    ki_matrix,    kind=wp)**2
      eps_i = cmplx(nr_inclusion, ki_inclusion, kind=wp)**2
      eps_eff = maxwell_garnett(eps_m, eps_i, f_inclusion)
      m_eff   = sqrt(eps_eff)
      if (aimag(m_eff) < 0.0_wp) m_eff = -m_eff
      nr_eff = real(m_eff, kind=wp)
      ki_eff = aimag(m_eff)
   end subroutine maxwell_garnett_m

end module effective_medium
