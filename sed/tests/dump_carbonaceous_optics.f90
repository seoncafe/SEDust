! Writes the carbonaceous grain optics of qpah (qpah_dl07 and qpah_ld01 for
! both charges, qpah_sca) on a grid of (a, lambda) that straddles the PAH
! cutoff (1/lambda = 17.25 um^-1), for each graphite source of the xi blend
! ('d16_sphere', 'd03_sphere', 'd16_spheroid') and for the two carbon-count
! settings (nc_coeff = 417 real-valued, the astrodust one; 470 rounded, the
! DL07 one), so that two library builds can be compared bit for bit.
! Compiled with -DCOMBINED_OPTICS it also writes what carbonaceous_grain_optics
! returns for both vintages at the same points.
!
! Not part of the library.  Build against a library directory LIB:
!   ifort -O2 -qopenmp -fpp [-DCOMBINED_OPTICS] -ILIB dump_carbonaceous_optics.f90 \
!         -LLIB -lsedust -L<hdf5>/lib -lhdf5_fortran -lhdf5 -lsz -lz -ldl -lm
! Usage: dump_carbonaceous_optics.x <data_dir> <output file>
program dump_carbonaceous_optics
   use sed_paths, only: sed_set_data_root
   use qpah,      only: qpah_dl07, qpah_ld01, qpah_sca, qpah_graphite_source, &
                        nc_coeff, nc_integer
#ifdef COMBINED_OPTICS
   use qpah,      only: carbonaceous_grain_optics
#endif
   implicit none
   integer, parameter :: dp = kind(1d0)
   integer, parameter :: NRAD = 24, NLAM_LOG = 120
   character(len=12), parameter :: SOURCES(3) = [character(len=12) :: &
        'd16_sphere', 'd03_sphere', 'd16_spheroid']
   character(len=512) :: data_dir, out_file
   real(dp), allocatable :: rad(:), lam(:)
   real(dp) :: q(10)
   real(dp) :: lam_cut
   integer :: u, isrc, inc, ja, jw, nlam

   call get_command_argument(1, data_dir)
   call get_command_argument(2, out_file)
   call sed_set_data_root(trim(data_dir))

   ! Radii 3.5e-4 .. 1 um (the PAH sizes and past the D16 sphere floor 1e-3 um).
   allocate(rad(NRAD))
   do ja = 1, NRAD
      rad(ja) = 3.5d-4 * (1.0d0/3.5d-4)**(real(ja-1, dp)/real(NRAD-1, dp))
   end do
   ! Wavelengths 1e-4 .. 1e3 um, log spaced, plus the cutoff itself and its
   ! two neighboring doubles.
   nlam = NLAM_LOG + 3
   allocate(lam(nlam))
   do jw = 1, NLAM_LOG
      lam(jw) = 1.0d-4 * (1.0d7)**(real(jw-1, dp)/real(NLAM_LOG-1, dp))
   end do
   lam_cut = 1.0d0/17.25d0
   lam(NLAM_LOG+1) = lam_cut
   lam(NLAM_LOG+2) = nearest(lam_cut,  1.0d0)
   lam(NLAM_LOG+3) = nearest(lam_cut, -1.0d0)

   open(newunit=u, file=trim(out_file), access='stream', status='replace')
   do isrc = 1, size(SOURCES)
      qpah_graphite_source = SOURCES(isrc)
      do inc = 1, 2
         if (inc == 1) then
            nc_coeff = 417.0d0;  nc_integer = .false.
         else
            nc_coeff = 470.0d0;  nc_integer = .true.
         end if
         do ja = 1, NRAD
            do jw = 1, nlam
               q = 0.0d0
               call qpah_dl07(0, rad(ja), lam(jw), q(1))
               call qpah_dl07(1, rad(ja), lam(jw), q(2))
               call qpah_ld01(0, rad(ja), lam(jw), q(3))
               call qpah_ld01(1, rad(ja), lam(jw), q(4))
               call qpah_sca(rad(ja), lam(jw), q(5), q(6))
               write(u) isrc, inc, ja, jw, rad(ja), lam(jw), q(1:6)
#ifdef COMBINED_OPTICS
               call carbonaceous_grain_optics('dl07', rad(ja), lam(jw), q(1), q(2), q(5), q(6))
               call carbonaceous_grain_optics('ld01', rad(ja), lam(jw), q(3), q(4), q(7), q(8))
               write(u) q(1:8)
#endif
            end do
         end do
      end do
   end do
   close(u)
end program dump_carbonaceous_optics
