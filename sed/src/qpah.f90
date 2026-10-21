module qpah
   ! PAH absorption cross sections per Draine & Li (2007), with the
   ! graphite transition active.
   !
   ! HD23 §3.2 / DL07 eq. 5-7 (= HD23 eq. 15-16):
   !   Q^PAH,Z(a, λ) = [1 - ξ_gra(a)] · Q^Z_PAH(a, λ) + ξ_gra(a) · Q^gra(a, λ)
   ! with Z = neutral or cation, and
   !   ξ_gra(a) = 0.01                                  for a ≤ 50 Å
   !   ξ_gra(a) = 0.01 + 0.99·[1 − (50 Å / a)³]          for a > 50 Å
   ! Below the DL07 PAH cutoff (λ ≤ 1/17.25 μm ≈ 0.058 μm), C^PAH is
   ! not defined and we use pure graphite (no ξ scaling), per
   ! the DL07 prescription.
   !
   ! Implements the QPAH_DL07 + drude_DL07 + cutoff
   ! trio; the graphite of the blend comes from whichever of three optics
   ! qpah_graphite_source names (see below), EXCEPT above the PAH cutoff,
   ! where the D16 tables' smallest radius (1e-3 um for the sphere table)
   ! is larger than the grains that matter and graphite is taken from the
   ! D03 dielectric functions in all three cases (see
   ! graphite_qabs_in_blend).
   ! Empirically (sed_total_d16_threeway.pdf)
   ! the D16 sphere agrees with HD23 PAH_irem.dat best in the
   ! 300-3000 um sub-mm band (PAH-only band-integrated ratio 0.99, against
   ! 1.05 for the D16 b/a=1.4 oblate spheroid and 0.61 for the D03 sphere),
   ! consistent with HD23 inheriting DL07's spherical-PAH framework.

   use constants,                 only: wp, pi
   use q_graphite_d16_sphere_mod, only: q_graphite_d16_sphere_abs
   use q_graphite_d16_mod,        only: q_graphite_d16_spheroid_abs => q_graphite_d16_abs
   use q_graphite_mod,            only: q_graphite_d03_abs => q_graphite_abs, &
                                       q_graphite_full
   implicit none
   private
   public :: qpah_dl07
   public :: qpah_ld01     ! Li & Draine 2001 PAH variant
   ! Which published PAH absorption cross section the carbonaceous blend takes:
   !   'dl07'  Draine & Li (2007), qpah_dl07 -- the default, and what the DL07
   !           model and the astrodust PAH population are defined with.  It
   !           carries the PAH cation near-infrared additions of DL07 eq. 2
   !           (the 1.05 and 1.26 um Drude resonances, the 1.905 um negative
   !           term and the NIR continuum, after Mattioda et al. 2005).
   !   'ld01'  Li & Draine (2001), qpah_ld01 -- the earlier vintage, with none
   !           of those, and with its own N_C = 468 (a/10A)^3 rather than
   !           nc_coeff.  Draine's 2003 kext_albedo table was computed with it.
   ! qpah_abs dispatches on this, so a caller that wants whichever vintage the
   ! model is configured for asks for it once rather than branching itself.
   public :: qpah_xsec_vintage
   public :: qpah_abs
   character(len=8), save :: qpah_xsec_vintage = 'dl07'
   ! Graphite optics used for the carbonaceous PAH<->graphite blend:
   !   'd16_sphere'    Draine 2016 turbostratic graphite, Maxwell-Garnett
   !                   effective medium, on a sphere (q_D16graphite.dat).
   !                   The HD23 choice and the production default.
   !   'd03_sphere'    Draine 2003 graphite dielectric functions with Mie on
   !                   a sphere, (1/3)E||c + (2/3)E_|_c. The original DL07
   !                   choice; what the dl07 model of calc_sed and calc_kext selects.
   !   'd16_spheroid'  the same D16 turbostratic material on Draine's
   !                   b/a = 1.4 oblate spheroid (qlib_gra_D16MGemt_1.400),
   !                   random-orientation averaged. Sensitivity test only.
   public :: qpah_graphite_source
   public :: xi_graphite, qpah_sca
   ! Absorption of both charge states and the scattering of one carbonaceous
   ! grain at one (a, lambda), with the D03 graphite Mie solved once for all
   ! three.  See carbonaceous_grain_optics.
   public :: carbonaceous_grain_optics
   ! Absorption of both charge states alone, with the graphite of the blend
   ! taken once for the pair -- or from a D03 graphite Q_abs the caller has
   ! already solved at the same (a, lambda).  See carbonaceous_absorption.
   public :: carbonaceous_absorption
   character(len=16), save :: qpah_graphite_source = 'd16_sphere'
   ! Carbon-atom-count coefficient Nc = nc_coeff*(a/10A)^3. Default 417 =
   ! HD23 eq.21 (rho_PAH = 2.0 g/cm^3), used by the astrodust path. Draine's
   ! DL07 prescription uses 470 (rho ~ 2.2; makeqlib NCAR = NINT(4.70e11*a_um^3))
   ! and rounds Nc to an integer -- the dl07 model sets nc_coeff=470, nc_integer=.true.
   ! to match Draine exactly on the DL07-model path.
   public :: nc_coeff, nc_integer
   real(kind=wp), save :: nc_coeff   = 417.0d0
   logical,       save :: nc_integer = .false.
   ! Toggle for empirical mode-by-mode sigma rescaling tuned to HD23
   ! PAH_irem.dat. Default .false. preserves the DL07 Table 1 values
   ! (with the Seon sigma_Ion(4) typo correction already in place).
   ! When .true., the multiplicative factors sigma_tune_neu(30) and
   ! sigma_tune_ion(30) below are applied to sigma_neu(j) and
   ! sigma_ion(j) inside drude_dl07.
   public :: use_tuned_pah, sigma_tune_neu, sigma_tune_ion, gamma_tune
   ! HD23 eq.15 graphite-blend parameters. Production defaults match
   ! DL07/HD23 prescription. Overridable at runtime for sensitivity
   ! scans.
   public :: xi_A_T_um, xi_FGMIN

   integer, parameter :: NMODE_PUB = 30
   logical, save :: use_tuned_pah = .false.
   ! Multiplicative correction factor for each mode; default 1.0 (no change).
   ! Set by calc_sed.f90 or by an external tuner before sed_solve.
   real(kind=wp), save :: sigma_tune_neu(NMODE_PUB) = 1.0_wp
   real(kind=wp), save :: sigma_tune_ion(NMODE_PUB) = 1.0_wp
   ! Multiplier for the Drude width gamma_j of each mode (DL07 Table 1
   ! column 3). Active only when use_tuned_pah=.true. Default 1.0.
   real(kind=wp), save :: gamma_tune(NMODE_PUB)    = 1.0_wp
   ! xi_gra(a) blending parameters (HD23 eq.16):
   !   xi_gra(a) = FGMIN                                    a <= A_T
   !             = FGMIN + (1-FGMIN) * [1 - (A_T/a)^3]      a >  A_T
   ! Defaults are production / HD23: A_T = 50 A = 5e-3 um, FGMIN = 0.01.
   real(kind=wp), save :: xi_A_T_um = 5.0e-3_wp
   real(kind=wp), save :: xi_FGMIN  = 0.01_wp

contains

   pure function xi_graphite(radius) result(xi)
      ! Graphite fraction of the carbonaceous cross section, HD23 eq. 16
      ! (= DL07 eq. 7): the transition from a PAH molecule to a graphite
      ! grain, complete well above A_T = 50 A and floored at FGMIN = 0.01
      ! below it.  Every cross section that carries the transition -- the
      ! two absorption vintages and the scattering -- reads it from here, so
      ! there is one curve and not three copies of it.
      real(kind=wp), intent(in) :: radius      ! [um]
      real(kind=wp) :: xi
      if (radius <= xi_A_T_um) then
         xi = xi_FGMIN
      else
         xi = xi_FGMIN + (1.0_wp - xi_FGMIN) * (1.0_wp - (xi_A_T_um/radius)**3)
      end if
   end function xi_graphite


   subroutine qpah_sca(radius, lambda, Qsca, gsca)
      ! Scattering of the carbonaceous grain that qpah_abs gives the
      ! absorption of.  A PAH molecule has no scattering prescription --
      ! DL07 and HD23 write C^neu and C^ion for absorption alone -- so what
      ! scatters is the graphite fraction of HD23 eq. 15,
      !
      !     Q_sca(a, lambda) = xi_gra(a) * Q_sca^gra(a, lambda) ,
      !
      ! and g is the graphite value, since graphite is the only thing in the
      ! mixture that scatters at all.  Above the PAH cutoff (21.4 eV) the
      ! carbonaceous grain IS graphite, exactly as it is in absorption there,
      ! so xi drops out.
      !
      ! The graphite comes from the D03 dielectric functions by Mie, the same
      ! random-orientation average (1/3 E||c + 2/3 E-perp-c) qpah_abs takes
      ! its 'd03_sphere' branch from.  The two D16 sources are absorption
      ! tables and carry no Q_sca or g, so they cannot serve this; the
      ! difference is a material one and it lives in the FIR, while the
      ! scattering this routine returns matters only in the far-UV.
      real(kind=wp), intent(in)  :: radius        ! [um]
      real(kind=wp), intent(in)  :: lambda        ! [um]
      real(kind=wp), intent(out) :: Qsca, gsca
      real(kind=wp) :: Qext_gra, Qsca_gra, Qabs_gra

      call q_graphite_full(radius, lambda, Qext_gra, Qsca_gra, Qabs_gra, gsca)
      Qsca = carbonaceous_qsca(radius, lambda, Qsca_gra)
   end subroutine qpah_sca


   pure function carbonaceous_qsca(radius, lambda, Qsca_gra) result(Qsca)
      ! HD23 eq. 15 for scattering: the graphite fraction xi_gra(a) scatters,
      ! the PAH part does not; above the PAH cutoff (1/lambda >= 17.25 um^-1)
      ! the grain is graphite and xi drops out.
      real(kind=wp), intent(in) :: radius        ! [um]
      real(kind=wp), intent(in) :: lambda        ! [um]
      real(kind=wp), intent(in) :: Qsca_gra      ! graphite (D03 Mie) Q_sca
      real(kind=wp) :: Qsca, x

      x = 1.0_wp / lambda
      if (x >= 17.25d0) then
         Qsca = Qsca_gra
      else
         Qsca = xi_graphite(radius) * Qsca_gra
      end if
   end function carbonaceous_qsca


   subroutine carbonaceous_grain_optics(vintage, radius, lambda, Qabs_neu, Qabs_ion, Qsca, gsca)
      ! The optics of one carbonaceous grain at one (a, lambda): the absorption
      ! efficiency of the neutral and of the cation (qpah_dl07 or qpah_ld01,
      ! as vintage names) and the scattering efficiency and asymmetry
      ! (qpah_sca).  Each output is bit for bit what that routine returns on
      ! its own; the gain is that the graphite Mie on the D03 dielectric
      ! functions (two polarizations) is solved ONCE here, where the three
      ! separate calls solve it up to three times -- once for the scattering
      ! at every wavelength, and once more for each charge state wherever the
      ! absorption takes D03 graphite (shortward of the PAH cutoff, 0.058 um,
      ! for every qpah_graphite_source; at all wavelengths for 'd03_sphere').
      ! That Mie is what the model build spends its time on.  The D16 table
      ! lookup, where the source is a D16 table, is likewise shared by the two
      ! charge states.
      character(len=*), intent(in)  :: vintage     ! 'dl07' | 'ld01'
      real(kind=wp),    intent(in)  :: radius      ! [um]
      real(kind=wp),    intent(in)  :: lambda      ! [um]
      real(kind=wp),    intent(out) :: Qabs_neu, Qabs_ion, Qsca, gsca
      real(kind=wp) :: Qext_d03, Qsca_d03, Qabs_d03

      call q_graphite_full(radius, lambda, Qext_d03, Qsca_d03, Qabs_d03, gsca)
      Qsca = carbonaceous_qsca(radius, lambda, Qsca_d03)
      call carbonaceous_absorption(vintage, radius, lambda, Qabs_neu, Qabs_ion, Qabs_d03)
   end subroutine carbonaceous_grain_optics


   subroutine carbonaceous_absorption(vintage, radius, lambda, Qabs_neu, Qabs_ion, Qabs_d03)
      ! Absorption efficiency of the neutral and of the cation carbonaceous
      ! grain at one (a, lambda), DL07 eq. 5-7 / HD23 eq. 15-16 with the PAH
      ! molecule of the vintage named (qpah_dl07 or qpah_ld01).  Each output is
      ! bit for bit what that routine returns on its own.  The graphite of the
      ! blend (graphite_qabs_in_blend) is the same for both charge states, so
      ! it is taken once here rather than once for each.
      !
      ! Qabs_d03, when present, is the D03 graphite Q_abs at this (a, lambda)
      ! that the caller has already solved (q_graphite_full: the same
      ! (1/3) E||c + (2/3) E-perp-c average of the same two Mie solutions
      ! q_graphite_d03_abs returns), and it is used wherever the blend takes
      ! D03 graphite -- shortward of the PAH cutoff (1/lambda >= 17.25 um^-1)
      ! for every qpah_graphite_source, and at all wavelengths for
      ! 'd03_sphere'.  A caller that tabulates the graphite and the
      ! carbonaceous grain at the same cells then solves that Mie once.
      character(len=*), intent(in)  :: vintage     ! 'dl07' | 'ld01'
      real(kind=wp),    intent(in)  :: radius      ! [um]
      real(kind=wp),    intent(in)  :: lambda      ! [um]
      real(kind=wp),    intent(out) :: Qabs_neu, Qabs_ion
      real(kind=wp),    intent(in), optional :: Qabs_d03
      real(kind=wp) :: Qabs_gra, x

      x = 1.0_wp / lambda
      if (present(Qabs_d03) .and. &
          (x >= 17.25d0 .or. trim(qpah_graphite_source) == 'd03_sphere')) then
         Qabs_gra = Qabs_d03
      else
         call graphite_qabs_in_blend(radius, lambda, Qabs_gra)
      end if

      select case (trim(vintage))
      case ('dl07')
         Qabs_neu = carbonaceous_qabs(radius, lambda, pah_molecule_qabs_dl07(0, radius, lambda), Qabs_gra)
         Qabs_ion = carbonaceous_qabs(radius, lambda, pah_molecule_qabs_dl07(1, radius, lambda), Qabs_gra)
      case ('ld01')
         Qabs_neu = carbonaceous_qabs(radius, lambda, pah_molecule_qabs_ld01(0, radius, lambda), Qabs_gra)
         Qabs_ion = carbonaceous_qabs(radius, lambda, pah_molecule_qabs_ld01(1, radius, lambda), Qabs_gra)
      case default
         write(*,'(a)') 'carbonaceous_absorption: unknown vintage "'//trim(vintage)//'"'
         write(*,'(a)') '      expected dl07 or ld01'
         stop 1
      end select
   end subroutine carbonaceous_absorption


   subroutine q_graphite_xi_blend_abs(radius, lambda, Qabs_gra)
      ! Graphite Q_abs entering the xi_gra blend of DL07 eq. 5-7 /
      ! HD23 eq. 15-16, taken from the optics qpah_graphite_source names.
      ! Every one of the three is a random-orientation-averaged absorption
      ! efficiency, which is what an unaligned carbonaceous population needs.
      ! All three underlying modules load their table on first use, and every
      ! caller reaches this from serial code, so no explicit load is issued
      ! here.
      real(kind=wp), intent(in)  :: radius        ! [um]
      real(kind=wp), intent(in)  :: lambda        ! [um]
      real(kind=wp), intent(out) :: Qabs_gra

      select case (trim(qpah_graphite_source))
      case ('d16_sphere')
         call q_graphite_d16_sphere_abs(radius, lambda, Qabs_gra)
      case ('d03_sphere')
         call q_graphite_d03_abs(radius, lambda, Qabs_gra)
      case ('d16_spheroid')
         call q_graphite_d16_spheroid_abs(radius, lambda, Qabs_gra)
      case default
         write(*,'(a)') 'qpah: unknown qpah_graphite_source "'// &
                        trim(qpah_graphite_source)//'"'
         write(*,'(a)') '      expected d16_sphere, d03_sphere or d16_spheroid'
         stop 1
      end select
   end subroutine q_graphite_xi_blend_abs


   subroutine graphite_qabs_in_blend(radius, lambda, Qabs_gra)
      ! Graphite Q_abs that DL07 eq. 5-7 / HD23 eq. 15-16 blend with the PAH
      ! one, for either vintage.  Above the PAH cutoff (21.4 eV) C^PAH = 0, so
      ! graphite IS the whole carbonaceous absorption -- and there it is
      ! needed down to a ~ 4e-4 um, below the D16 sphere table's smallest
      ! radius (D16_LAM_MIN_UM's companion floor a = 1e-3 um).  That table
      ! clamps to its boundary radius, and since Q_abs is in the Rayleigh
      ! limit (x = 2 pi a / lambda ~ 0.04) where Q_abs ~ a, the clamp
      ! overestimates the absorption by (1e-3 um)/a, up to 2.4x.  So graphite
      ! is taken from the D03 dielectric functions there: Mie on the
      ! random-orientation average (1/3 E||c + 2/3 E_|_c) with no radius
      ! clamp, valid to 6.2e-5 um.  This holds for every setting of
      ! qpah_graphite_source, since it is the tables' radius floor and not a
      ! preference that forces it.  Longward of the cutoff the selected
      ! source takes over: there xi_gra weights graphite by only 0.01 for the
      ! small grains, and the D16 sphere is what matches HD23 PAH_irem.dat in
      ! the sub-mm (see the module header).
      ! The D03 branch is unreachable unless the caller's wavelength grid is
      ! carried below 1/17.25 um = 0.05797 um.  The T-matrix Q table stops at
      ! 0.0912 um (1/0.0912 = 10.96 < 17.25), so a run on the plain table grid
      ! never enters it.
      real(kind=wp), intent(in)  :: radius        ! [um]
      real(kind=wp), intent(in)  :: lambda        ! [um]
      real(kind=wp), intent(out) :: Qabs_gra
      real(kind=wp) :: x

      x = 1.0_wp / lambda
      if (x >= 17.25d0) then
         call q_graphite_d03_abs(radius, lambda, Qabs_gra)   ! DL07 D03 graphite
      else
         call q_graphite_xi_blend_abs(radius, lambda, Qabs_gra)
      end if
   end subroutine graphite_qabs_in_blend


   pure function carbonaceous_qabs(radius, lambda, Qabs_pah, Qabs_gra) result(Qabs)
      ! DL07 eq. 5-7 / HD23 eq. 15-16: the carbonaceous Q_abs as the
      ! xi_gra(a)-weighted blend of the PAH molecule and graphite.  Shortward
      ! of the PAH cutoff (1/lambda >= 17.25 um^-1) no PAH C_abs is defined
      ! and the grain is pure graphite.
      real(kind=wp), intent(in) :: radius        ! [um]
      real(kind=wp), intent(in) :: lambda        ! [um]
      real(kind=wp), intent(in) :: Qabs_pah, Qabs_gra
      real(kind=wp) :: Qabs, x, xi_gra

      x = 1.0_wp / lambda
      xi_gra = xi_graphite(radius)
      if (x >= 17.25d0) then
         Qabs = Qabs_gra
      else
         Qabs = (1.0_wp - xi_gra) * Qabs_pah + xi_gra * Qabs_gra
      end if
   end function carbonaceous_qabs


   subroutine qpah_abs(charge, radius, lambda, Qabs)
      ! Carbonaceous absorption efficiency in whichever published vintage
      ! qpah_xsec_vintage names.  An unrecognized value stops the run rather
      ! than silently selecting one of the two.
      integer,       intent(in)  :: charge        ! 0 = neutral, 1 = cation
      real(kind=wp), intent(in)  :: radius        ! [um]
      real(kind=wp), intent(in)  :: lambda        ! [um]
      real(kind=wp), intent(out) :: Qabs
      select case (trim(qpah_xsec_vintage))
      case ('dl07');  call qpah_dl07(charge, radius, lambda, Qabs)
      case ('ld01');  call qpah_ld01(charge, radius, lambda, Qabs)
      case default
         write(*,'(a)') 'qpah: unknown qpah_xsec_vintage "'// &
                        trim(qpah_xsec_vintage)//'"'
         write(*,'(a)') '      expected dl07 or ld01'
         stop 1
      end select
   end subroutine qpah_abs


   subroutine qpah_dl07(charge, radius, lambda, Qabs)
      ! Draine & Li (2007) carbonaceous absorption efficiency: the PAH
      ! molecule blended with graphite by xi_gra(a).
      integer,       intent(in)  :: charge        ! 0 = neutral, 1 = cation
      real(kind=wp), intent(in)  :: radius        ! [um]
      real(kind=wp), intent(in)  :: lambda        ! [um]
      real(kind=wp), intent(out) :: Qabs
      real(kind=wp) :: Qabs_gra

      call graphite_qabs_in_blend(radius, lambda, Qabs_gra)
      Qabs = carbonaceous_qabs(radius, lambda, pah_molecule_qabs_dl07(charge, radius, lambda), Qabs_gra)
   end subroutine qpah_dl07


   function pah_molecule_qabs_dl07(charge, radius, lambda) result(Qabs_pah)
      ! Absorption efficiency of the PAH molecule alone, C^PAH N_C / (pi a^2),
      ! with C^PAH the DL07 cross section for one carbon atom (the continuum
      ! pieces below plus the Drude terms of DL07 Table 1).
      ! Zero shortward of the PAH cutoff, where C^PAH is not defined.
      integer,       intent(in)  :: charge        ! 0 = neutral, 1 = cation
      real(kind=wp), intent(in)  :: radius        ! [um]
      real(kind=wp), intent(in)  :: lambda        ! [um]
      real(kind=wp) :: Qabs_pah

      real(kind=wp) :: CPAH
      real(kind=wp) :: x, H_C
      real(kind=wp) :: Nc
      integer :: jmode

      ! Number of carbon atoms, Nc = nc_coeff*(a/10A)^3. nc_coeff defaults to
      ! 417 = HD23 eq.21 (rho_PAH = 2.0 g/cm^3) for the astrodust path; the
      ! DL07 path sets nc_coeff=470 (rho ~ 2.2) and nc_integer=.true.
      ! Since Qabs = C^PAH * Nc / (pi a^2), Qabs scales LINEARLY with Nc, so
      ! the 417->470 difference is a uniform ~12.7% on the carbonaceous Qabs.
      Nc = nc_coeff * (radius/1.0d-3)**3
      if (nc_integer) Nc = anint(Nc)

      ! H/C ratio — DL07 eq. (4.7).
      if (Nc <= 25d0) then
         H_C = 0.5d0
      else if (Nc < 100d0) then
         H_C = 0.5d0 * sqrt(25d0/Nc)
      else
         H_C = 0.25d0
      end if

      x = 1d0 / lambda                            ! lambda^-1 [1/um]

      if (x >= 17.25d0) then
         CPAH = 0.0d0
      else if (x >= 15.0d0) then
         CPAH = (126.0d0 - 6.4943d0*x) * 1d-18
      else if (x >= 10.0d0) then
         CPAH = drude_dl07(charge, 1, lambda, H_C) + (-3.0d0 + 1.35d0*x) * 1d-18
      else if (x >= 7.7d0) then
         CPAH = (66.302d0 + x*(-24.367d0 + x*(2.95d0 - 0.1057d0*x))) * 1d-18
      else if (x >= 5.9d0) then
         CPAH = drude_dl07(charge, 2, lambda, H_C) + &
                (1.8687d0 + 0.1905d0*x + 0.4175d0*(x-5.9d0)**2 + &
                 0.04370d0*(x-5.9d0)**3) * 1.d-18
      else if (x >= 3.3d0) then
         CPAH = drude_dl07(charge, 2, lambda, H_C) + &
                (1.8687d0 + 0.1905d0*x) * 1.d-18
      else
         ! x < 3.3 (lambda > 0.303 um) — IR / mid-IR features
         CPAH = 34.58d0 * 10d0**(-18.0d0 - 3.431d0/x) * cutoff(charge, Nc, lambda)
         do jmode = 3, 30
            CPAH = CPAH + drude_dl07(charge, jmode, lambda, H_C)
         end do
      end if

      ! Cation near-IR continuum, DL07 eq. 2 (Mattioda et al. 2005a):
      !   C_abs/N_C = 3.5e-19 * 10^(-1.45/x) * exp[-(0.1x)^2] cm^2, x = 1/lambda[um].
      ! The Gaussian factor carries 0.1x inside the square, as in Draine.
      if (charge /= 0 .and. x < 17.25d0) then
         CPAH = CPAH + 3.5d0 * 10d0**(-19d0 - 1.45d0/x) * exp(-(0.1d0*x)**2)
      end if

      ! C^PAH is the cross section per C atom [cm^2]; convert to Q via Q = C·Nc/(πa²).
      Qabs_pah = CPAH * Nc / (pi * (radius*1d-4)**2)
   end function pah_molecule_qabs_dl07


   function drude_dl07(charge, i, lambda, H_C) result(D)
      ! Drude profile contribution at mode index i (DL07 Table 1).
      integer,       intent(in)  :: charge, i
      real(kind=wp), intent(in)  :: lambda, H_C
      real(kind=wp) :: D
      integer, parameter :: nmode = 30
      real(kind=wp), parameter :: lambdaj(nmode) = [ &
         0.0722d0, 0.2175d0, 1.050d0,  1.260d0,  1.905d0,  3.300d0,  5.270d0,  5.700d0,  6.220d0,  6.690d0, &
         7.417d0,  7.598d0,  7.850d0,  8.330d0,  8.610d0, 10.68d0,  11.23d0,  11.33d0,  11.99d0,  12.62d0, &
        12.69d0,  13.48d0,  14.19d0,  15.90d0,  16.45d0,  17.04d0,  17.375d0, 17.87d0,  18.92d0,  15.0d0   ]
      real(kind=wp), parameter :: width(nmode) = [ &
         0.195d0,  0.217d0,  0.055d0,  0.11d0,   0.09d0,   0.012d0,  0.034d0,  0.035d0,  0.030d0,  0.070d0, &
         0.126d0,  0.044d0,  0.053d0,  0.052d0,  0.039d0,  0.020d0,  0.012d0,  0.032d0,  0.045d0,  0.042d0, &
         0.013d0,  0.040d0,  0.025d0,  0.020d0,  0.014d0,  0.065d0,  0.012d0,  0.016d0,  0.10d0,   0.8d0    ]
      real(kind=wp), parameter :: sigma_neu(nmode) = [ &
         7.97d7,   1.23d7,   0.0d0,    0.0d0,    0.0d0,  394.0d0,    2.5d0,    4.0d0,   29.4d0,    7.35d0, &
        20.8d0,   18.1d0,   21.9d0,    6.94d0,  27.8d0,   0.3d0,   18.9d0,   52.0d0,   24.2d0,   35.0d0,   &
         1.3d0,    8.0d0,    0.45d0,   0.04d0,   0.5d0,   2.22d0,   0.11d0,   0.067d0,  0.10d0,  50.0d0    ]
      real(kind=wp), parameter :: sigma_ion(nmode) = [ &
         7.97d7,   1.23d7,   2.0d4,    7.8d3, -146.5d0,  89.4d0,   20.0d0,   32.0d0,  235.0d0,   59.0d0,   &
       181.0d0,  163.0d0,  197.0d0,   48.0d0, 194.0d0,   0.3d0,   17.7d0,   49.0d0,   20.5d0,   31.0d0,    &
         1.3d0,    8.0d0,    0.45d0,   0.04d0,   0.5d0,   2.22d0,   0.11d0,   0.067d0,  0.17d0,  50.0d0    ]
      real(kind=wp) :: ratio, sigma

      ! 1d-4 converts um to cm.
      if (charge == 0) then
         sigma = sigma_neu(i) * 1d-20 * 1d-4
         if (use_tuned_pah) sigma = sigma * sigma_tune_neu(i)
      else
         sigma = sigma_ion(i) * 1d-20 * 1d-4
         if (use_tuned_pah) sigma = sigma * sigma_tune_ion(i)
      end if

      ! H/C scaling on H-related modes (3.3, 11.3, 11.9, ..., 14.2 um group)
      if (i == 6 .or. (i >= 14 .and. i <= 22)) sigma = sigma * H_C

      ratio = lambda / lambdaj(i)
      if (use_tuned_pah) then
         block
            real(kind=wp) :: w
            w = width(i) * gamma_tune(i)
            D = 2d0/pi * (lambdaj(i)*w*sigma) / ((ratio - 1d0/ratio)**2 + w**2)
         end block
      else
         D = 2d0/pi * (lambdaj(i)*width(i)*sigma) / ((ratio - 1d0/ratio)**2 + width(i)**2)
      end if
   end function drude_dl07


   function cutoff(charge, Nc, lambda) result(c)
      ! Desert et al. (1990) cutoff function.
      integer,       intent(in) :: charge
      real(kind=wp), intent(in) :: Nc, lambda
      real(kind=wp) :: c, M, lambda_c, y

      ! M = number of fused benzenoid rings.
      if (Nc >= 40d0) then
         M = 0.4d0 * Nc
      else
         M = 0.3d0 * Nc
      end if

      ! Cutoff wavelength (um), DL07 eq. 4.5
      if (charge == 0) then
         lambda_c = 1.0d0 / (3.804d0/sqrt(M) + 1.052d0)
      else
         lambda_c = 1.0d0 / (2.282d0/sqrt(M) + 0.889d0)
      end if
      y = lambda_c / lambda
      c = atan(1000.0d0 * (y - 1.0d0)**3 / y) / pi + 0.5d0
   end function cutoff

   ! ------------------------------------------------------------------
   ! Li & Draine (2001) carbonaceous absorption (QPAH_LD01 ported verbatim
   ! of DL07). The IR/UV continuum pieces are identical to
   ! qpah_dl07; the differences are (i) the Drude profiles use the LD01
   ! Table-1 sigmas with the E6.2/E7.7/E8.6 enhancements and H/C scaling
   ! (drude_ld01, 14 modes), (ii) Nc = 468 (a/nm)^3, and (iii) NO DL07
   ! cation continuum boost. Blended with graphite by the same xi_gra(a).
   ! ------------------------------------------------------------------
   subroutine qpah_ld01(charge, radius, lambda, Qabs)
      integer,       intent(in)  :: charge        ! 0 = neutral, 1 = cation
      real(kind=wp), intent(in)  :: radius        ! [um]
      real(kind=wp), intent(in)  :: lambda        ! [um]
      real(kind=wp), intent(out) :: Qabs
      real(kind=wp) :: Qabs_gra

      call graphite_qabs_in_blend(radius, lambda, Qabs_gra)
      Qabs = carbonaceous_qabs(radius, lambda, pah_molecule_qabs_ld01(charge, radius, lambda), Qabs_gra)
   end subroutine qpah_ld01


   function pah_molecule_qabs_ld01(charge, radius, lambda) result(Qabs_pah)
      ! Absorption efficiency of the PAH molecule alone in the Li & Draine
      ! (2001) vintage, C^PAH N_C / (pi a^2) with N_C = 468 (a/10 A)^3.
      integer,       intent(in)  :: charge        ! 0 = neutral, 1 = cation
      real(kind=wp), intent(in)  :: radius        ! [um]
      real(kind=wp), intent(in)  :: lambda        ! [um]
      real(kind=wp) :: Qabs_pah
      real(kind=wp) :: CPAH, x, H_C, Nc
      integer :: jmode

      Nc = 468.0d0 * (radius/1.0d-3)**3
      if (Nc <= 25d0) then
         H_C = 0.5d0
      else if (Nc < 100d0) then
         H_C = 0.5d0 / sqrt(Nc/25d0)
      else
         H_C = 0.25d0
      end if

      x = 1d0 / lambda
      if (x >= 17.25d0) then
         CPAH = 0.0d0
      else if (x >= 15.0d0) then
         CPAH = (126.0d0 - 6.4943d0*x) * 1d-18
      else if (x >= 10.0d0) then
         CPAH = drude_ld01(charge, 1, lambda, H_C) + (-3.0d0 + 1.35d0*x) * 1d-18
      else if (x >= 7.7d0) then
         CPAH = (66.302d0 + x*(-24.367d0 + x*(2.95d0 - 0.1057d0*x))) * 1d-18
      else if (x >= 5.9d0) then
         CPAH = drude_ld01(charge, 2, lambda, H_C) + &
                (1.8687d0 + 0.1905d0*x + 0.4175d0*(x-5.9d0)**2 + &
                 0.04370d0*(x-5.9d0)**3) * 1.d-18
      else if (x >= 3.3d0) then
         CPAH = drude_ld01(charge, 2, lambda, H_C) + &
                (1.8687d0 + 0.1905d0*x) * 1.d-18
      else
         CPAH = 34.58d0 * 10d0**(-18.0d0 - 3.431d0/x) * cutoff(charge, Nc, lambda)
         do jmode = 3, 14
            CPAH = CPAH + drude_ld01(charge, jmode, lambda, H_C)
         end do
      end if

      Qabs_pah = CPAH * Nc / (pi * (radius*1d-4)**2)
   end function pah_molecule_qabs_ld01

   function drude_ld01(charge, i, lambda, H_C) result(D)
      ! LD01 Table 1 Drude profile for mode i (sigma in 1e-20 cm^2/C).
      integer,       intent(in) :: charge, i
      real(kind=wp), intent(in) :: lambda, H_C
      real(kind=wp) :: D
      real(kind=wp), parameter :: lambdaj(14) = &
         (/0.0722d0,0.2175d0,3.3d0,6.2d0,7.7d0,8.6d0,11.3d0,11.9d0,12.7d0, &
           16.4d0,18.3d0,21.2d0,23.1d0,26.0d0/)
      real(kind=wp), parameter :: width(14) = &
         (/0.195d0,0.217d0,0.012d0,0.030d0,0.091d0,0.047d0,0.018d0,0.025d0, &
           0.024d0,0.010d0,0.036d0,0.038d0,0.046d0,0.69d0/)
      real(kind=wp), parameter :: sigma_Neu(14) = &
         (/7.97d7,1.23d7,197.d0,19.6d0,60.9d0,34.7d0,427.d0,72.7d0, &
           167.d0,5.52d0,6.04d0,10.8d0,2.78d0,15.2d0/)
      real(kind=wp), parameter :: sigma_Ion(14) = &
         (/7.97d7,1.23d7,44.7d0,157.d0,548.d0,242.d0,400.d0,61.4d0, &
           149.d0,5.52d0,6.04d0,10.8d0,2.78d0,15.2d0/)
      real(kind=wp), parameter :: E6_2 = 3.0d0, E7_7 = 2.0d0, E8_6 = 2.0d0
      real(kind=wp) :: ratio, sigma

      if (charge == 0) then
         sigma = sigma_Neu(i)*1d-20*1d-4   ! 1e-4: um -> cm
      else
         sigma = sigma_Ion(i)*1d-20*1d-4
      end if
      if (i == 4) sigma = sigma * E6_2
      if (i == 5) sigma = sigma * E7_7
      if (i == 6) sigma = sigma * E8_6
      if (i == 3 .or. i == 6) sigma = sigma * H_C
      if (i == 7 .or. i == 8 .or. i == 9) sigma = sigma * H_C/3d0
      ratio = lambda / lambdaj(i)
      D = 2d0/pi * (lambdaj(i)*width(i)*sigma) / ((ratio - 1d0/ratio)**2 + width(i)**2)
   end function drude_ld01

end module qpah
