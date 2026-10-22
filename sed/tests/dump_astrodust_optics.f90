! Writes the astrodust model's grain optics, as build_dust leaves them, to one
! binary file, so that two library builds can be compared bit for bit.
!
! The build is the one an RT host makes for a photoionization run (MoCHII's
! grain_model_mod): astrodust, NT = 200, T = 2.7 .. 5000 K, the EUV axis of the
! HDF5 product, sd_index = 7, u_isrf = 1.  Written, in this order, each as
! (name, n1, n2, values):
!   module arrays of sed_astrodust_mod  Cabs_cneu, Cabs_cion, Cabs_pah,
!                                       kappB_pah_first, kappCMB_pah
!   every population of the model       Cabs, Csca, gsca, kappB, kappCMB, dn
!   size_integrated_extinction          Cext, Cabs, Csca, gbar, albedo
! and the wall time of build_dust on standard output.
!
! Not part of the library.  Build against a library directory LIB:
!   ifort -O2 -qopenmp -fpp -ILIB dump_astrodust_optics.f90 -LLIB -lsedust \
!         -L<hdf5>/lib -lhdf5_fortran -lhdf5 -lsz -lz -ldl
! Usage: dump_astrodust_optics.x <data_dir> <output file> [ld01]
! A third argument 'ld01' sets use_ld01_pah_xsec, the Li & Draine (2001) PAH
! cross sections of the two PAH populations.
program dump_astrodust_optics
   use dust_lib
   use sed_astrodust_mod, only: Cabs_cneu, Cabs_cion, Cabs_pah, kappB_pah_first, &
                                kappCMB_pah
   use omp_lib,           only: omp_get_wtime
   use pah_ld01_mod,      only: use_ld01_pah_xsec
   implicit none
   integer, parameter :: dp = kind(1d0)
   type(dust_model_t) :: m
   character(len=512) :: data_dir, out_file, msg, vintage
   character(len=16)  :: tag
   integer :: st, u, ip, nlam
   real(dp) :: t0, t_build
   real(dp), allocatable :: Cext(:), Cabs(:), Csca(:), gbar(:), albedo(:)

   call get_command_argument(1, data_dir)
   call get_command_argument(2, out_file)
   call get_command_argument(3, vintage)
   use_ld01_pah_xsec = trim(vintage) == 'ld01'

   t0 = omp_get_wtime()
   call build_dust(m, 'astrodust', trim(data_dir), 200, 2.7d0, 5.0d3, include_euv=.true., &
        status=st, message=msg, sd_index=7, u_isrf=1.0d0)
   t_build = omp_get_wtime() - t0
   if (st /= 0) then
      write(*,*) 'build_dust status=', st, ' : ', trim(msg);  stop 2
   end if
   write(*,'(a,f10.3)') 'build_dust wall time [s] = ', t_build

   nlam = m%NLAM
   allocate(Cext(nlam), Cabs(nlam), Csca(nlam), gbar(nlam), albedo(nlam))
   call size_integrated_extinction(m, Cext, Cabs, Csca, gbar=gbar, albedo=albedo, status=st)
   if (st /= 0) then
      write(*,*) 'size_integrated_extinction status=', st;  stop 3
   end if

   open(newunit=u, file=trim(out_file), access='stream', status='replace')
   call put2('Cabs_cneu', Cabs_cneu)
   call put2('Cabs_cion', Cabs_cion)
   call put2('Cabs_pah', Cabs_pah)
   call put2('kappB_pah_first', kappB_pah_first)
   call put1('kappCMB_pah', kappCMB_pah)
   do ip = 1, size(m%pops)
      write(tag, '(a,i0,a)') 'pop', ip, '_'
      call put2(trim(tag)//'Cabs', m%pops(ip)%Cabs)
      call put2(trim(tag)//'Csca', m%pops(ip)%Csca)
      call put2(trim(tag)//'gsca', m%pops(ip)%gsca)
      call put2(trim(tag)//'kappB', m%pops(ip)%kappB)
      call put1(trim(tag)//'kappCMB', m%pops(ip)%kappCMB)
      call put1(trim(tag)//'dn', m%pops(ip)%dn)
   end do
   call put1('ext_Cext', Cext)
   call put1('ext_Cabs', Cabs)
   call put1('ext_Csca', Csca)
   call put1('ext_gbar', gbar)
   call put1('ext_albedo', albedo)
   close(u)

contains

   subroutine put2(name, a)
      character(len=*), intent(in) :: name
      real(dp),         intent(in) :: a(:,:)
      character(len=16) :: nm
      nm = name
      write(u) nm, size(a, 1), size(a, 2), a
   end subroutine put2

   subroutine put1(name, a)
      character(len=*), intent(in) :: name
      real(dp),         intent(in) :: a(:)
      character(len=16) :: nm
      nm = name
      write(u) nm, size(a), 1, a
   end subroutine put1

end program dump_astrodust_optics
