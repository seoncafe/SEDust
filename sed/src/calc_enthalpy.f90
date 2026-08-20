program calc_enthalpy
   ! Sanity / debug tool: tabulate U(T, a_eff) for the two astrodust
   ! enthalpy stages and write them as ASCII tables.
   !
   ! NOT required to compute an SED. The SED driver (calc_sed.f90) calls
   ! the module functions enthalpy_S1 / enthalpy_S2 directly in its
   ! inner loop; closed-form Debye sums are microseconds per call and
   ! a precomputed table just adds I/O + grid-mismatch risk.
   !
   ! Use this tool to:
   !   - compare Stage 2 vs Stage 1 enthalpy at fixed (T, a)
   !   - plot U(T) curves for documentation
   !   - spot-check against published DL01 reference values
   !
   ! Stages:
   !   S1 - Stage 1, silicate-only, literal DL01 prefactor (rho=3.5)
   !   S2 - Stage 2, silicate + carbonaceous, volume-weighted
   !
   ! Grids:
   !   T  - 201 log-spaced points from 2.7 K (CMB floor) to 5000 K
   !   a  - 169 points read from data/astrodust/DH21_aeff (matches Q grid)
   !
   ! Output (run from sed/, files dropped in sed/output/):
   !   output/enthalpy_S1.dat
   !   output/enthalpy_S2.dat
   ! Each: 4 header lines + (NT * NA) data rows of "T[K]  a_eff[um]  U[erg]",
   ! T outer / a inner, mirroring the q_table layout.
   !
   ! Usage:
   !   ./calc_enthalpy.x [c2]
   !
   ! `c2` is the one axis this program has a referent for: the Stage-1
   ! density-corrected prefactor of enthalpy_S1.  It tags its table
   ! (output/enthalpy_S1_c2.dat) so the default one is not overwritten, and no
   ! Stage-2 table is written for it -- enthalpy_S2 does not read the setting,
   ! so a tagged copy would be a second name for the same numbers.  Every other
   ! axis -- solver, radiation field, wavelength grid -- is refused by name: an
   ! enthalpy table has no referent for one.

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   use enthalpy_astrodust_mod, only: enthalpy_S1, enthalpy_S2, s1_density_corrected
   use sed_run_options, only: run_options_t, declare_run_options, &
                              read_run_option, check_run_options, &
                              run_options_tag, write_run_option_usage
   implicit none

   character(len=*), parameter :: F_AEFF  = '../data/astrodust/DH21_aeff'

   integer,  parameter :: NA = 169, NT = 201
   real(wp), parameter :: T_LO = 2.7_wp, T_HI = 5.0e3_wp

   real(wp) :: a_eff(NA), T_grid(NT)
   integer  :: i, ja, jt, narg, iarg
   real(wp) :: U
   integer  :: u_s1, u_s2
   logical  :: taken, write_s2
   type(run_options_t) :: o
   character(len=64)   :: arg
   character(len=160)  :: tag
   character(len=64)   :: f_out_s1, f_out_s2

   ! ---- the one axis this program has a referent for ---------------------
   call declare_run_options(o, program='calc_enthalpy', stage1=.true.)
   narg = command_argument_count()
   do iarg = 1, narg
      call get_command_argument(iarg, arg)
      call read_run_option(trim(arg), o, taken)
      if (.not. taken) then
         write(*,'(a,a)') ' calc_enthalpy: unknown argument ', trim(arg)
         write(*,'(a)')   ' usage: ./calc_enthalpy.x [c2]'
         call write_run_option_usage(o)
         stop 1
      end if
   end do
   call check_run_options(o)
   ! This program's one setting lands in the module it already links, so it
   ! needs no applier of its own.
   s1_density_corrected = o%stage1_density_corrected
   tag = run_options_tag(o)
   ! Stage 2 does not read the Stage-1 prefactor, so a tagged Stage-2 table
   ! would be a second name for the numbers the default one already carries.
   write_s2 = (len_trim(tag) == 0)
   if (len_trim(tag) > 0) then
      f_out_s1 = 'output/enthalpy_S1_'//trim(tag)//'.dat'
      f_out_s2 = ''
   else
      f_out_s1 = 'output/enthalpy_S1.dat'
      f_out_s2 = 'output/enthalpy_S2.dat'
   end if

   ! Read a_eff grid
   call read_one_col(F_AEFF, NA, a_eff)
   write(*,'(a,i0,a)') ' a_eff: ', NA, ' values'
   write(*,'(a,2es12.4,a)') '   range = ', a_eff(1), a_eff(NA), ' [um]'

   ! Build T grid
   do i = 1, NT
      T_grid(i) = T_LO * (T_HI / T_LO)**(real(i-1, wp) / real(NT-1, wp))
   end do
   write(*,'(a,i0,a)') ' T:     ', NT, ' values (log-spaced)'
   write(*,'(a,2es12.4,a)') '   range = ', T_grid(1), T_grid(NT), ' [K]'

   ! S1
   open(newunit=u_s1, file=trim(f_out_s1), status='replace')
   write(u_s1,'(a)') '# DH21 astrodust enthalpy table'
   write(u_s1,'(a)') '# Stage 1 (silicate-only, literal DL01 prefactor, rho=3.5)'
   write(u_s1,'(a,i0,a,i0)') '# NT = ', NT, '   NA = ', NA
   write(u_s1,'(a)') '# T[K]            a_eff[um]       U[erg]'

   if (write_s2) then
      open(newunit=u_s2, file=trim(f_out_s2), status='replace')
      write(u_s2,'(a)') '# DH21 astrodust enthalpy table'
      write(u_s2,'(a)') '# Stage 2 (silicate + carbonaceous, volume-weighted)'
      write(u_s2,'(a,i0,a,i0)') '# NT = ', NT, '   NA = ', NA
      write(u_s2,'(a)') '# T[K]            a_eff[um]       U[erg]'
   end if

   do jt = 1, NT
      do ja = 1, NA
         U = enthalpy_S1(T_grid(jt), a_eff(ja))
         write(u_s1,'(es16.8,1x,es12.4,1x,es16.8)') T_grid(jt), a_eff(ja), U
         if (write_s2) then
            U = enthalpy_S2(T_grid(jt), a_eff(ja))
            write(u_s2,'(es16.8,1x,es12.4,1x,es16.8)') T_grid(jt), a_eff(ja), U
         end if
      end do
   end do
   close(u_s1)
   if (write_s2) close(u_s2)

   write(*,'(a,a)') ' wrote ', trim(f_out_s1)
   if (write_s2) then
      write(*,'(a,a)') ' wrote ', trim(f_out_s2)
   else
      write(*,'(a)') ' Stage 2 does not read this setting; no Stage-2 table written.'
   end if

contains

   subroutine read_one_col(filename, n, x)
      ! DH21_aeff: 2 header lines + 1 long data line.
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: n
      real(wp),         intent(out) :: x(n)
      integer :: u, ios
      character(len=512) :: header
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a,a)') 'ERROR: cannot open ', trim(filename); stop 1
      end if
      read(u,'(a)') header
      read(u,'(a)') header
      read(u,*) x(1:n)
      close(u)
   end subroutine read_one_col

end program calc_enthalpy
