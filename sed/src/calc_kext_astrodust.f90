program calc_kext_astrodust
   !====================================================================
   ! Size-distribution-integrated extinction cross section per H, albedo,
   ! and scattering asymmetry <cos> for the HD23 astrodust+PAH Milky-Way
   ! model, computed from EXACTLY the optics the SED pipeline uses.
   !
   ! We call sed_init() -- the same initializer main_astrodust.x uses --
   ! which loads the random-orientation T-matrix Q table
   ! (q_astrodust_P0.20_Fe0.00_1.400.dat) and the HD23 size distribution,
   ! and builds the single-grain cross sections:
   !   astrodust grains : Cabs(NLAM,NA), Csca(NLAM,NA) from the T-matrix
   !                      Q table (Qabs, Qsca x pi a^2); the asymmetry g
   !                      is taken from the same table (q_table_mod::gpar)
   !                      and interpolated onto the size grid the same way.
   !   PAH              : charge-resolved DL07 absorption Cabs_cneu/Cabs_cion
   !                      (neutral + cation), mixed by f_ion with the 100 A
   !                      cutoff exactly as sed_init builds dn_cneu/dn_cion.
   !                      The pipeline does not model PAH scattering
   !                      (negligible Rayleigh), so in this model the PAH
   !                      contributes to extinction through absorption only
   !                      and the astrodust grains carry all the scattering
   !                      and the asymmetry.
   !
   ! Size integrals per H atom (dn from size_distribution.dat are already binned):
   !   C_abs/H = sum_a [ dn_Ad Cabs_Ad + dn_cneu Cabs_neu + dn_cion Cabs_ion ]
   !   C_sca/H = sum_a   dn_Ad Csca_Ad
   !   C_ext/H = C_abs/H + C_sca/H
   !   albedo  = C_sca / C_ext
   !   <cos>   = ( sum_a dn_Ad Csca_Ad g_Ad ) / C_sca
   !
   ! Output: ../data/kext_astrodust_MW.dat  (kext_albedo-style columns).
   !====================================================================
   use constants,         only: wp, pi
   use sed_astrodust_mod, only: build_astrodust, dust_extinction, dust_model_t, &
                                NLAM, NA, lam, aeff, dn_ad, dn_pah
   use enthalpy_astrodust_mod, only: RHO_AD     ! astrodust bulk density [g/cm^3] = 2.74
   implicit none

   ! PAH mass density (HD23 eq. 21 / nc_coeff=417 convention).
   real(wp), parameter :: RHO_PAH = 2.0_wp      ! [g/cm^3]
   real(wp) :: Mdust_H, acm3                     ! dust mass per H [g/H]

   character(len=*), parameter :: F_QT  = '../tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: F_SD  = '../data/release/size_distribution.dat'
   character(len=*), parameter :: F_EXT = '../data/release/extinction.dat'
   character(len=*), parameter :: F_SCA = '../data/release/scattering.dat'
   character(len=*), parameter :: F_OUT = '../data/kext_astrodust_MW.dat'

   type(dust_model_t)    :: m
   real(wp), allocatable :: Cext(:), Cabs_t(:), Csca_t(:), alb(:), gbar(:)
   integer  :: jw, ja, u

   write(*,'(a)') ' calc_kext_astrodust: building the astrodust model ...'
   ! Build the model the RT hosts use, so this table and their in-memory optics
   ! come from one calculation.  (NT/T_lo/T_hi only affect the thermal tables,
   ! which the extinction does not use.)
   call build_astrodust(m, F_QT, F_SD, 100, 1.0_wp, 3000.0_wp)
   write(*,'(a,i0,a,i0,a)') '   NLAM=', NLAM, '  NA=', NA, '  (size grid)'

   ! Dust mass per H [g/H], self-consistent with the size distribution and
   ! the model grain densities: M = (4/3) pi a^3 rho, summed over the binned
   ! number per H of astrodust (rho_Ad=2.74) and PAH (rho_PAH=2.0) grains.
   Mdust_H = 0.0_wp
   do ja = 1, NA
      acm3 = (aeff(ja) * 1.0e-4_wp)**3
      Mdust_H = Mdust_H + (4.0_wp/3.0_wp)*pi*acm3 * &
                ( dn_ad(ja)*RHO_AD + dn_pah(ja)*RHO_PAH )
   end do
   write(*,'(a,es12.5,a)') '   M_dust/H = ', Mdust_H, ' g/H'

   ! Size integral over every population of the model: the astrodust grains
   ! carry the scattering and its asymmetry, the PAHs enter through absorption.
   allocate(Cext(NLAM), Cabs_t(NLAM), Csca_t(NLAM), alb(NLAM), gbar(NLAM))
   call dust_extinction(m, Cext, Cabs_t, Csca_t, gbar=gbar, albedo=alb)

   ! ---- write the table -------------------------------------------------
   open(newunit=u, file=F_OUT, status='replace', action='write')
   write(u,'(a)') '# Extinction, albedo, and scattering asymmetry for the'
   write(u,'(a)') '# Hensley & Draine (2023) astrodust+PAH Milky-Way model.'
   write(u,'(a)') '#'
   write(u,'(a)') '# Astrodust grains: random-orientation T-matrix optics'
   write(u,'(a)') '#   (Draine & Hensley 2021 dielectric, P=0.20, fFe=0, b/a=1.4);'
   write(u,'(a)') '#   Cabs, Csca and <cos> from the T-matrix Q table.'
   write(u,'(a)') '# PAH: charge-resolved DL07 absorption (neutral+cation by f_ion,'
   write(u,'(a)') '#   100 A -> 100% ionized); PAH scattering is negligible and not'
   write(u,'(a)') '#   modelled, so PAH adds to extinction via absorption only.'
   write(u,'(a)') '# Size integral over the HD23 release size distribution (per H).'
   write(u,'(a)') '# Computed from the SED pipeline optics (sed_init), all components.'
   write(u,'(a,es13.6,a)') '# Dust mass per H, M_dust/N_H = ', Mdust_H, &
        ' g/H  (rho_Ad=2.74, rho_PAH=2.0 g/cm^3; K_abs = C_abs/H / this).'
   write(u,'(a)') '#'
   write(u,'(a)') '#   lambda      albedo      <cos>      C_ext/H        C_abs/H' // &
                  '        C_sca/H         K_abs'
   write(u,'(a)') '#  (micron)                            (cm^2/H)       (cm^2/H)' // &
                  '       (cm^2/H)       (cm^2/g)'
   do jw = 1, NLAM
      write(u,'(es13.5e3,2(1x,f10.6),4(1x,es15.7e3))') &
         lam(jw), alb(jw), gbar(jw), Cext(jw), Cabs_t(jw), Csca_t(jw), &
         Cabs_t(jw)/Mdust_H
   end do
   close(u)
   write(*,'(a,a)') ' wrote ', F_OUT

   ! ---- validation against the HD23 release (total = Ad + PAH) ----------
   call validate()

contains

   subroutine validate()
      ! Compare C_ext/H and C_sca/H against the HD23 release total columns
      ! (extinction.dat / scattering.dat: lambda, tau_Ad, tau_PAH, tau_tot).
      real(wp), allocatable :: wr(:), ext_tot(:), sca_tot(:)
      real(wp) :: a(4)
      integer  :: ur, ios, n, i
      character(len=512) :: line
      real(wp) :: ce, cs, dext, dsca
      real(wp), parameter :: bands(5) = (/0.15_wp, 0.55_wp, 2.2_wp, 12.0_wp, 100.0_wp/)
      n = 0
      open(newunit=ur, file=F_EXT, status='old', action='read', iostat=ios)
      if (ios /= 0) then; write(*,'(a)') ' (release not found; skipping validation)'; return; end if
      do
         read(ur,'(a)',iostat=ios) line; if (ios /= 0) exit
         line = adjustl(line); if (len_trim(line)==0 .or. line(1:1)=='#') cycle
         n = n + 1
      end do
      rewind(ur)
      allocate(wr(n), ext_tot(n), sca_tot(n))
      i = 0
      do
         read(ur,'(a)',iostat=ios) line; if (ios /= 0) exit
         line = adjustl(line); if (len_trim(line)==0 .or. line(1:1)=='#') cycle
         i = i + 1; read(line,*) a; wr(i) = a(1); ext_tot(i) = a(4)
      end do
      close(ur)
      open(newunit=ur, file=F_SCA, status='old', action='read', iostat=ios)
      i = 0
      do
         read(ur,'(a)',iostat=ios) line; if (ios /= 0) exit
         line = adjustl(line); if (len_trim(line)==0 .or. line(1:1)=='#') cycle
         i = i + 1; read(line,*) a; sca_tot(i) = a(4)
      end do
      close(ur)
      write(*,'(a)') ' validation vs HD23 release (total Ad+PAH):'
      write(*,'(a)') '   lam[um]   C_ext ours/HD23   C_sca ours/HD23'
      do i = 1, 5
         ce = loginterp(lam, Cext,   NLAM, bands(i))
         cs = loginterp(lam, Csca_t, NLAM, bands(i))
         dext = loginterp(wr, ext_tot, n, bands(i))
         dsca = loginterp(wr, sca_tot, n, bands(i))
         write(*,'(f10.3,2(7x,f10.4))') bands(i), ce/dext, cs/dsca
      end do
   end subroutine validate

   function loginterp(x, y, n, x0) result(y0)
      integer,  intent(in) :: n
      real(wp), intent(in) :: x(n), y(n), x0
      real(wp) :: y0
      integer  :: j
      do j = 2, n
         if (x(j) >= x0) exit
      end do
      if (j > n) j = n
      y0 = exp( log(max(y(j-1),tiny(0.0_wp))) + &
           (log(x0)-log(x(j-1)))/(log(x(j))-log(x(j-1))) * &
           (log(max(y(j),tiny(0.0_wp)))-log(max(y(j-1),tiny(0.0_wp)))) )
   end function loginterp

end program calc_kext_astrodust
