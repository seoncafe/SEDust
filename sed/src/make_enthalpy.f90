program make_enthalpy
   ! Sanity / debug tool: tabulate U(T, a_eff) for the two astrodust
   ! enthalpy stages and write them as ASCII tables.
   !
   ! NOT required to compute an SED. The SED driver (main_astrodust.f90) calls
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
   !   a  - 169 points read from data/dielectric/DH21_aeff (matches Q grid)
   !
   ! Output (run from sed/, files dropped in sed/output/):
   !   output/enthalpy_S1.dat
   !   output/enthalpy_S2.dat
   ! Each: 4 header lines + (NT * NA) data rows of "T[K]  a_eff[um]  U[erg]",
   ! T outer / a inner, mirroring the q_table layout.
   !
   ! Usage:
   !   ./make_enthalpy.x

   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   use enthalpy_astrodust_mod, only: enthalpy_S1, enthalpy_S2
   implicit none

   character(len=*), parameter :: F_AEFF  = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: F_OUT_S1 = 'output/enthalpy_S1.dat'
   character(len=*), parameter :: F_OUT_S2 = 'output/enthalpy_S2.dat'

   integer,  parameter :: NA = 169, NT = 201
   real(wp), parameter :: T_LO = 2.7_wp, T_HI = 5.0e3_wp

   real(wp) :: a_eff(NA), T_grid(NT)
   integer  :: i, ja, jt
   real(wp) :: U
   integer  :: u_s1, u_s2

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
   open(newunit=u_s1, file=F_OUT_S1, status='replace')
   write(u_s1,'(a)') '# DH21 astrodust enthalpy table'
   write(u_s1,'(a)') '# Stage 1 (silicate-only, literal DL01 prefactor, rho=3.5)'
   write(u_s1,'(a,i0,a,i0)') '# NT = ', NT, '   NA = ', NA
   write(u_s1,'(a)') '# T[K]            a_eff[um]       U[erg]'

   open(newunit=u_s2, file=F_OUT_S2, status='replace')
   write(u_s2,'(a)') '# DH21 astrodust enthalpy table'
   write(u_s2,'(a)') '# Stage 2 (silicate + carbonaceous, volume-weighted)'
   write(u_s2,'(a,i0,a,i0)') '# NT = ', NT, '   NA = ', NA
   write(u_s2,'(a)') '# T[K]            a_eff[um]       U[erg]'

   do jt = 1, NT
      do ja = 1, NA
         U = enthalpy_S1(T_grid(jt), a_eff(ja))
         write(u_s1,'(es16.8,1x,es12.4,1x,es16.8)') T_grid(jt), a_eff(ja), U
         U = enthalpy_S2(T_grid(jt), a_eff(ja))
         write(u_s2,'(es16.8,1x,es12.4,1x,es16.8)') T_grid(jt), a_eff(ja), U
      end do
   end do
   close(u_s1); close(u_s2)

   write(*,'(a,a)') ' wrote ', F_OUT_S1
   write(*,'(a,a)') ' wrote ', F_OUT_S2

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

end program make_enthalpy
