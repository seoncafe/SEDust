program main_dl07
   ! Draine & Li (2007) model driver: amorphous silicate + carbonaceous
   ! (PAH + graphite) grains with WD01 size distributions, heated by the
   ! Mathis ISRF at U = 1. Reproduces the DL07spec reference spectra.
   !
   ! Usage:  ./main_dl07.x [model] [euv] [qm|draine]
   !   model = mw31_00 | mw31_10 | ... | mw31_60 (default) | mw31_avg-style,
   !           lmc2_00 | lmc2_05 | lmc2_10 | smc
   ! Output: output/dl07_sed_ours_<model>[_euv][_<solver>].dat
   !   columns: lambda[um]  lamI_total/NH  lamI_sil/NH  lamI_carb/NH
   !            [erg s^-1 cm^-2 sr^-1 H^-1] (the HD23/DL07 convention)
   !   `euv` solves on the _euv Q table, whose grid runs to 1.0e-4 um instead
   !   of stopping at the Lyman limit; the Mathis field carries no photon there.

   use, intrinsic :: iso_fortran_env, only: real64
   use constants,         only: wp
   use radfield,          only: J_Mathis, use_mathis_corrected
   use sed_astrodust_mod, only: sed_init_dl07, sed_solve_dl07, NLAM, lam, stoch_method
   use qpah,              only: qpah_use_d03_graphite, nc_coeff, nc_integer
   use grain_dist_mod,    only: gd_apply_d03_reduction
   implicit none

   ! DL07 / WD MW R_V=3.1_60 = the 2003 model: WD01 abundances reduced by the
   ! Draine-2003a factor 0.93 (gd_apply_d03_reduction=.true.), heated by the
   ! MMP83 field. The astrodust "mathis_corrected" field (w_4000=1.65e-13) IS
   ! the canonical MMP83 (matches u_star=8.64e-13 to 0.4%); the uncorrected
   ! literal 1e-13 under-normalizes the optical band by ~7%. So use the
   ! corrected field. (Earlier 0.93-off + uncorrected was an accidental
   ! cancellation that desynced emission from the 2003 extinction.)

   ! The DL07 optics are Mie on the D03 dielectric functions, but the model's
   ! wavelength grid still comes from the astrodust Q table, so which of the two
   ! shipped tables this driver reads decides how far the SED reaches.  Default
   ! is the non-ionizing one: the Mathis field is zero below the Lyman limit, so
   ! the 633 extra wavelengths of the _euv table carry no photon.
   character(len=*), parameter :: F_QTAB =  &
      '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: F_QTAB_EUV =  &
      '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400_euv.dat'
   character(len=*), parameter :: F_SIZE = '../data/release/size_distribution.dat'
   real(wp),         parameter :: U_ISRF = 1.0_wp        ! DL07 reference: U = 1
   real(wp),         parameter :: T_LO = 2.7_wp, T_HI = 5.0e3_wp
   integer,          parameter :: NT_IN = 200
   character(len=*), parameter :: OUTDIR = 'output/dl07_sed_ours_'

   real(wp), allocatable :: J_lam(:)
   real(wp), allocatable :: lamI_tot(:), lamI_sil(:), lamI_carb(:)
   integer  :: sd_index, k, u, narg, iarg
   logical  :: use_qm, use_draine
   ! Solve over the ionizing band as well, on the _euv Q table.
   logical  :: use_euv_grid = .false.
   character(len=32) :: model, arg
   character(len=64) :: fname
   character(len=8)  :: qmtag
   character(len=8)  :: euvtag

   ! ---- Args (any order): model name, and optional 'qm' to use the
   ! energy-space transition-matrix (thermal-discrete) stochastic solver
   ! instead of the default GD (Draine-narrowing) solver. ----
   model  = 'mw31_60'
   use_qm = .false.
   use_draine = .false.
   narg = command_argument_count()
   do iarg = 1, narg
      call get_command_argument(iarg, arg)
      if (trim(arg) == 'qm') then
         use_qm = .true.
      else if (trim(arg) == 'draine') then
         use_draine = .true.        ! Draine's original GD (default is 'heuristic')
      else if (trim(arg) == 'euv') then
         use_euv_grid = .true.      ! solve on the _euv Q table's 1762 grid
      else
         model = trim(arg)
      end if
   end do

   select case (trim(model))
   case ('mw31_00'); sd_index = 1
   case ('mw31_10'); sd_index = 2
   case ('mw31_20'); sd_index = 3
   case ('mw31_30'); sd_index = 4
   case ('mw31_40'); sd_index = 5
   case ('mw31_50'); sd_index = 6
   case ('mw31_60'); sd_index = 7
   case ('lmc2_00'); sd_index = 29
   case ('lmc2_05'); sd_index = 30
   case ('lmc2_10'); sd_index = 31
   case ('smc');     sd_index = 32
   case default
      write(*,'(a,a)') ' main_dl07: unknown model ', trim(model)
      write(*,'(a)')   ' valid: mw31_00..mw31_60, lmc2_00, lmc2_05, lmc2_10, smc'
      stop 1
   end select

   write(*,'(a)') '=========================================================='
   write(*,'(a)') ' main_dl07: Draine & Li (2007) silicate + carbonaceous SED'
   write(*,'(a)') '=========================================================='
   write(*,'(a,a)')   ' model         : ', trim(model)
   write(*,'(a,i0)')  ' WD01 index    : ', sd_index
   write(*,'(a,f6.3)')' U (Mathis)    : ', U_ISRF
   use_mathis_corrected  = .true.    ! corrected = canonical MMP83 (matches u_star 0.4%)
   qpah_use_d03_graphite = .true.    ! DL07 carbonaceous uses D03 graphite
   gd_apply_d03_reduction = .true.   ! 0.93 abundance reduction = the 2003 model
   nc_coeff  = 470.0d0               ! DL07 Nc coefficient (rho~2.2)
   nc_integer = .true.               ! Nc is rounded to an integer
   qmtag = '';  euvtag = ''
   if (use_euv_grid) euvtag = '_euv'
   if (use_qm) then
      stoch_method = 'qm'            ! energy-space transition matrix (dbdis)
      qmtag = '_qm'
   else if (use_draine) then
      stoch_method = 'draine'        ! Draine's original GD (default is 'heuristic')
      qmtag = '_draine'
   end if
   write(*,'(a,l1)')    ' mathis_corr   : ', use_mathis_corrected
   write(*,'(a,l1)')    ' d03_graphite  : ', qpah_use_d03_graphite
   write(*,'(a,l1)')    ' d03_reduction : ', gd_apply_d03_reduction
   write(*,'(a,f6.1,a,l1)') ' Nc_coeff      : ', nc_coeff, '   Nc_integer: ', nc_integer
   write(*,'(a,a)')     ' solver        : ', trim(stoch_method)

   if (use_euv_grid) then
      call sed_init_dl07(F_QTAB_EUV, F_SIZE, sd_index, U_ISRF, NT_IN, T_LO, T_HI)
   else
      call sed_init_dl07(F_QTAB, F_SIZE, sd_index, U_ISRF, NT_IN, T_LO, T_HI)
   end if
   write(*,'(a,i0,a)') ' sed_init_dl07 done. NLAM=', NLAM, '.'

   allocate(J_lam(NLAM), lamI_tot(NLAM), lamI_sil(NLAM), lamI_carb(NLAM))
   call J_Mathis(U_ISRF, lam, J_lam)

   write(*,'(a)') ' solving DL07 SED (silicate + carbonaceous) ...'
   call sed_solve_dl07(J_lam, lamI_tot, lamI_sil, lamI_carb)

   write(fname,'(a,a,a,a,a)') OUTDIR, trim(model), trim(euvtag), trim(qmtag), '.dat'
   open(newunit=u, file=trim(fname), status='replace', action='write')
   write(u,'(a,a,a,f5.2)') '# DL07 model SED (this work), ', trim(model), &
        ', Mathis ISRF U = ', U_ISRF
   write(u,'(a)') '# silicate (D03) + carbonaceous (DL07 PAH Nc=470 + D03 graphite blend)'
   write(u,'(a)') '# PAH ionization computed via WD01b grain charging (pah_ionfrac)'
   write(u,'(a)') '# columns: lambda[um]  lamI_total/NH  lamI_sil/NH  lamI_carb/NH'
   write(u,'(a)') '#          [erg s^-1 cm^-2 sr^-1 H^-1]'
   ! e3 on the intensities: on a grid carried into the EUV, lambda*I_lambda
   ! underflows to subnormals, and a two-digit exponent field drops the E there
   ! (4.29970510-319), which readers outside Fortran cannot parse.
   do k = 1, NLAM
      write(u,'(es14.6,3(1x,es16.8e3))') lam(k), lamI_tot(k), lamI_sil(k), lamI_carb(k)
   end do
   close(u)
   write(*,'(a,a)') ' wrote ', trim(fname)
   write(*,'(a)') ' main_dl07: done.'

   deallocate(J_lam, lamI_tot, lamI_sil, lamI_carb)
end program main_dl07
