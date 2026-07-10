module enthalpy_astrodust_mod
   ! Astrodust enthalpy U(T, a_eff).
   !
   ! Two stages, both used end-to-end:
   !
   !   Stage 1 (S1) - silicate-only (literal DL01 'Sil' Debye)
   !     n_form(a) = 5.1040126e10 * a^3 [um^3]
   !
   !   Stage 2 (S2) - silicate + carbonaceous, volume-weighted
   !     Astrodust solid = (1-P) * V is split into two sub-volumes by
   !     volume fraction f_C^vol (carbonaceous) and 1 - f_C^vol (silicate).
   !     Each sub-volume contributes its own DL01 Debye U.
   !     U_Ad = U_Sil(N_Sil) + U_Car(N_Car), additive.
   !     HD23 sec 2.3.1 gives f_C^vol ~ 0.19.
   !
   ! Implementation note. enthalpy_DL01 takes (T, R) and
   ! internally computes natom = (constant) * R^3. To get sub-volume
   ! atom counts we feed it a *scaled* radius
   !   R_eff = R * (volume fraction)^(1/3)
   ! which is mathematically equivalent (modulo the (natom-2) -2 term,
   ! negligible for any realistic grain size).

   use constants,    only: wp
   use enthalpy,     only: enthalpy_DL01
   implicit none
   private
   public :: enthalpy_S1, enthalpy_S2
   public :: P_PORO, F_C_VOL, F_SIL_VOL, RHO_AD
   public :: s1_density_corrected

   ! Stage-1 prefactor option (paper Fig. 4 / Sect. 5.1):
   !  .false. = C1 (literal DL01 astrosilicate, rho = 3.5 g/cm^3);
   !  .true.  = C2 (density-corrected to the bulk astrodust rho = 2.74,
   !            i.e. N_atom scaled by 2.74/3.5 = 0.7829).
   logical, save :: s1_density_corrected = .false.

   ! Reference parameters
   real(wp), parameter :: P_PORO       = 0.20_wp     ! porosity
   real(wp), parameter :: RHO_AD       = 2.74_wp     ! astrodust solid mass density [g/cm^3]
                                                     ! = rho_solid * (1-P), HD23 sec 2.3.2

   ! Carbonaceous volume fraction of the astrodust solid.
   ! Derived from DH21a Table 2 (fFe = 0 column):
   !   V_car          = 0.86e-27 cm^3 H^-1
   !   V_Ad*(1-P)     = 4.46e-27 cm^3 H^-1   (V_sil + V_mix + V_car + V_Fe inc)
   !   f_C^vol        = V_car / V_Ad(1-P) = 0.86 / 4.46 = 0.1928
   ! Stage 2 lumps V_sil + V_mix + V_Fe inc into the "silicate-like"
   ! sub-volume, on the convention that "everything not carbonaceous is
   ! treated as silicate", giving:
   !   f_Sil^vol = 1 - f_C^vol = 0.8072
   ! (Cross-check vs DH21a eq. 8: V_sil/V_solid = 0.52 — that's just the
   ! pure silicate fraction, not the lumped silicate-like.)
   ! Insensitive to fFe: at fFe=0.10, V_car/V_Ad(1-P) = 0.86/4.40 = 0.1955.
   real(wp), parameter :: F_C_VOL   = 0.193_wp
   real(wp), parameter :: F_SIL_VOL = 1.0_wp - F_C_VOL

contains

   function enthalpy_S1(T, radius_um) result(U)
      ! Stage 1: silicate-only DL01. C1 (literal, rho=3.5) by default;
      ! C2 (density-corrected to rho_Ad=2.74) when s1_density_corrected.
      ! C2 scales N_atom by 2.74/3.5, achieved by a radius rescaling
      ! since N_atom propto a^3: a_eff = a * (2.74/3.5)^(1/3).
      real(wp), intent(in) :: T, radius_um
      real(wp) :: U, r_eff

      if (s1_density_corrected) then
         r_eff = radius_um * (RHO_AD / 3.5_wp)**(1.0_wp/3.0_wp)
      else
         r_eff = radius_um
      end if
      U = enthalpy_DL01(T, r_eff, 'Sil ')
   end function enthalpy_S1


   function enthalpy_S2(T, radius_um) result(U)
      ! Stage 2: silicate + carbonaceous, volume-weighted on the
      ! (1-P)-corrected solid sub-volumes.
      !   N_Sil = (1-P)*(1-f_C^vol) * (literal Sil prefactor) * a^3
      !   N_Car = (1-P)* f_C^vol    * (literal Car0 prefactor) * a^3
      ! achieved by passing scaled radii to DL01.
      real(wp), intent(in) :: T, radius_um
      real(wp) :: U, R_sil_eff, R_car_eff, U_sil, U_car

      R_sil_eff = radius_um * ((1.0_wp - P_PORO) * F_SIL_VOL)**(1.0_wp/3.0_wp)
      R_car_eff = radius_um * ((1.0_wp - P_PORO) * F_C_VOL  )**(1.0_wp/3.0_wp)
      U_sil     = enthalpy_DL01(T, R_sil_eff, 'Sil ')
      U_car     = enthalpy_DL01(T, R_car_eff, 'Car0')
      U         = U_sil + U_car
   end function enthalpy_S2

end module enthalpy_astrodust_mod
