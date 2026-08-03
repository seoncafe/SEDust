module fixed_orientation_phase_matrix
   !! Phase (Mueller) matrix of a particle in a fixed orientation, formed from
   !! its 2x2 complex amplitude matrix.  This is the angular counterpart of
   !! fixed_orientation_amplitude: AMPL returns the amplitude matrix, this
   !! routine turns it into the matrix that transforms the Stokes vector.
   !!
   !! The sixteen bilinear combinations are Mishchenko's Z11...Z44 of
   !! ampld.lp.f, kept uncompiled in reference/upstream/ -- Eqs. (13)-(29) of
   !! "Calculation of the amplitude matrix for a nonspherical particle in a
   !! fixed orientation", Appl. Opt. 39, 1026-1031, 2000.  There they are
   !! inline statements of the main program that prints the matrix; here they
   !! are a named routine that writes it into a caller-supplied array.  No
   !! expression, operation order, or working precision differs from the
   !! original.  The routine holds no state, needs no workspace, and may be
   !! called concurrently.
   !!
   !! STOKES CONVENTION.  S = [[S11,S12],[S21,S22]] = [[VV,VH],[HV,HH]] in the
   !! meridional basis of each propagation direction: e_1 = theta-hat (V),
   !! e_2 = phi-hat (H), real unit vectors, so that Q = I_v - I_h.  That is
   !! the basis AMPL rotates its amplitudes into, and Z is the phase matrix of
   !! the Stokes vector (I, Q, U, V) written in it.  The physical element is
   !! the real part of each bilinear, which is what the routine returns.
   !!
   !! UNITS.  Z carries the square of the units of S.  AMPL returns amplitudes
   !! with the dimension of length, so with its wavelength in microns Z is a
   !! differential scattering cross section in um^2 sr^-1.
   !!
   !! The grain orientation enters only through S: the caller fixes the Euler
   !! angles when it evaluates the amplitude matrix, so an orientation average
   !! accumulates Z from this routine at whatever (ALPHA, BETA) it needs.
   !!
   !! Re-exported by tmatrix_api next to AMPL, which supplies its argument.
   use tmatrix_kinds, only: wp
   implicit none
   private
   public :: phase_matrix_from_amplitude

contains

   subroutine phase_matrix_from_amplitude(s11, s12, s21, s22, z)
      ! INPUT
      !   s11, s12, s21, s22   amplitude matrix elements, S = [[VV,VH],[HV,HH]]
      ! OUTPUT
      !   z(4,4)               phase matrix, in the square of the units of S
      complex(wp), intent(in)  :: s11, s12, s21, s22
      real(wp),    intent(out) :: z(4,4)
      complex(wp), parameter :: CI = (0.0_wp, 1.0_wp)

      z(1,1) = 0.5_wp*real( s11*conjg(s11)+s12*conjg(s12) &
                           +s21*conjg(s21)+s22*conjg(s22), kind=wp)
      z(1,2) = 0.5_wp*real( s11*conjg(s11)-s12*conjg(s12) &
                           +s21*conjg(s21)-s22*conjg(s22), kind=wp)
      z(1,3) = real(-s11*conjg(s12)-s22*conjg(s21), kind=wp)
      z(1,4) = real(CI*(s11*conjg(s12)-s22*conjg(s21)), kind=wp)

      z(2,1) = 0.5_wp*real( s11*conjg(s11)+s12*conjg(s12) &
                           -s21*conjg(s21)-s22*conjg(s22), kind=wp)
      z(2,2) = 0.5_wp*real( s11*conjg(s11)-s12*conjg(s12) &
                           -s21*conjg(s21)+s22*conjg(s22), kind=wp)
      z(2,3) = real(-s11*conjg(s12)+s22*conjg(s21), kind=wp)
      z(2,4) = real(CI*(s11*conjg(s12)+s22*conjg(s21)), kind=wp)

      z(3,1) = real(-s11*conjg(s21)-s22*conjg(s12), kind=wp)
      z(3,2) = real(-s11*conjg(s21)+s22*conjg(s12), kind=wp)
      z(3,3) = real( s11*conjg(s22)+s12*conjg(s21), kind=wp)
      z(3,4) = real(-CI*(s11*conjg(s22)+s21*conjg(s12)), kind=wp)

      z(4,1) = real(CI*(s21*conjg(s11)+s22*conjg(s12)), kind=wp)
      z(4,2) = real(CI*(s21*conjg(s11)-s22*conjg(s12)), kind=wp)
      z(4,3) = real(-CI*(s22*conjg(s11)-s12*conjg(s21)), kind=wp)
      z(4,4) = real( s22*conjg(s11)-s12*conjg(s21), kind=wp)
   end subroutine phase_matrix_from_amplitude

end module fixed_orientation_phase_matrix
