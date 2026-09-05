program test_dustem_product
   ! Does a DustEM-defined model built from its HDF5 product equal the same
   ! model built from the DustEM text tables?
   !
   !   ./test_dustem_product.x [themis | g18d | all]     (default: all)
   !
   ! data/<model>/sedust_<model>.h5 does not carry a rounded copy of a
   ! calculation.  It carries the published Q_/G_ tables interpolated onto each
   ! population's model radii -- the nsize points the GRAIN line spaces evenly
   ! in ln a between its amin and amax -- which is precisely the step
   ! build_dustem performs at build time when it reads the text tables itself,
   ! by the same routine.  That interpolation is in radius alone, so cutting
   ! the wavelength axis and interpolating commute, so the two routes follow
   ! the same steps.  What separates them is storage: the product stores its
   ! computed quantities as 32-bit (see the storage-precision note in
   ! sedust_h5.f90), whose relative resolution is 1.2e-7, while the DustEM text
   ! tables the other route reads carry 7 to 13 significant digits.  The two
   ! cannot agree better than float32 resolution, so the tolerance here is 1e-6
   ! relative -- the smallest round number safely above it -- for C_ext, C_abs,
   ! C_sca and <cos theta>.  Anything larger than that IS a disagreement in the
   ! interpolation, which is what this checks.
   !
   ! Whether the model HAS an asymmetry parameter is part of the model, so it
   ! is compared too: the G18D product carries no g dataset for the two
   ! populations the DustEM distribution ships no G_ file for, and a product
   ! route that wrote a zero there instead would report an asymmetry this model
   ! does not have.
   !
   ! v1.00 makes the same comparison inside its check_build_dust.x, which
   ! already builds every model both ways; this version has no such program,
   ! so the comparison stands on its own here.
   use constants,         only: wp
   use sed_astrodust_mod, only: dust_model_t, build_dustem, &
                                size_integrated_extinction, dust_mass_per_H, &
                                sed_verbose
   implicit none

   character(len=*), parameter :: F_GR_TH = '../data/themis/GRAIN_J13.DAT'
   character(len=*), parameter :: D_TH    = '../data/themis/'
   character(len=*), parameter :: H5_TH   = '../data/themis/sedust_themis.h5'
   character(len=*), parameter :: F_GR_G1 = '../data/g18d/GRAIN_G17_ModelD.DAT'
   character(len=*), parameter :: D_G1    = '../data/g18d/'
   character(len=*), parameter :: H5_G1   = '../data/g18d/sedust_g18d.h5'
   integer,  parameter :: NT_IN = 100
   real(wp), parameter :: T_LO = 1.0_wp, T_HI = 3000.0_wp
   real(wp), parameter :: TOL = 1.0e-6_wp

   character(len=32) :: which
   integer :: narg, nbad

   sed_verbose = .false.
   narg = command_argument_count()
   which = 'all'
   if (narg >= 1) call get_command_argument(1, which)

   nbad = 0
   select case (trim(which))
   case ('themis');  call check_model('themis', F_GR_TH, D_TH, H5_TH)
   case ('g18d');    call check_model('g18d',   F_GR_G1, D_G1, H5_G1)
   case ('all');     call check_model('themis', F_GR_TH, D_TH, H5_TH)
                     call check_model('g18d',   F_GR_G1, D_G1, H5_G1)
   case default
      write(*,'(a)') ' usage: ./test_dustem_product.x [themis | g18d | all]'
      stop 1
   end select

   write(*,'(a)') ''
   if (nbad == 0) then
      write(*,'(a)') ' PASS'
   else
      write(*,'(i0,a)') nbad, ' CHECK(S) FAILED'
      stop 1
   end if

contains

   subroutine check_model(model, grain, mdir, h5)
      character(len=*), intent(in) :: model, grain, mdir, h5
      call compare_routes(model//', non-EUV', grain, mdir, h5, .false.)
      call compare_routes(model//', EUV',     grain, mdir, h5, .true.)
   end subroutine check_model


   subroutine compare_routes(label, grain, mdir, h5, wide)
      character(len=*), intent(in) :: label, grain, mdir, h5
      logical,          intent(in) :: wide
      type(dust_model_t) :: mh, mt
      real(wp), allocatable :: eh(:), ah(:), sh(:), gh(:)
      real(wp), allocatable :: et(:), at(:), st_(:), gt(:)
      real(wp) :: dlam, mdh, mdt
      character(len=256) :: gm_h, gm_t
      integer :: st, n, ip
      logical :: bad

      write(*,'(a)') ''
      write(*,'(a,a)') ' === ', label

      call build_dustem(mh, grain, mdir, NT_IN, T_LO, T_HI, status=st, &
                        include_euv=wide, gsca_missing=gm_h, qtable_path=h5)
      call fail_if(st /= 0, 'product route: build_dustem status', st)
      call build_dustem(mt, grain, mdir, NT_IN, T_LO, T_HI, status=st, &
                        include_euv=wide, gsca_missing=gm_t)
      call fail_if(st /= 0, 'text route: build_dustem status', st)
      if (mh%NLAM /= mt%NLAM) then
         write(*,'(a,i0,a,i0)') '   NLAM differs: product ', mh%NLAM, ', text ', mt%NLAM
         nbad = nbad + 1;  return
      end if
      n = mh%NLAM
      dlam = maxval(abs(mh%lam/mt%lam - 1.0_wp))
      write(*,'(a,i0,a,es9.2)') '   NLAM = ', n, ',  lambda max rel ', dlam
      call fail_if(dlam > TOL, 'lambda', 0)

      allocate(eh(n), ah(n), sh(n), gh(n), et(n), at(n), st_(n), gt(n))
      call size_integrated_extinction(mh, eh, ah, sh, gbar=gh)
      call size_integrated_extinction(mt, et, at, st_, gbar=gt)
      call report('C_ext ', eh, et)
      call report('C_abs ', ah, at)
      call report('C_sca ', sh, st_)
      call report('<cos> ', gh, gt)
      deallocate(eh, ah, sh, gh, et, at, st_, gt)

      mdh = dust_mass_per_H(mh);  mdt = dust_mass_per_H(mt)
      write(*,'(a,es14.7,a,es14.7,a,es9.2)') '   M_dust/H  product ', mdh, &
           '  text ', mdt, '   rel ', abs(mdh/mdt - 1.0_wp)
      call fail_if(abs(mdh/mdt - 1.0_wp) > TOL, 'M_dust/H', 0)

      write(*,'(a,l1,a,l1)') '   gsca_complete  product ', mh%gsca_complete, &
           '  text ', mt%gsca_complete
      if (len_trim(gm_t) > 0) write(*,'(a,a)') '     no g published for: ', trim(gm_t)
      call fail_if(mh%gsca_complete .neqv. mt%gsca_complete, 'gsca_complete', 0)
      call fail_if(trim(gm_h) /= trim(gm_t), 'which populations report no g', 0)
      bad = size(mh%pops) /= size(mt%pops)
      if (.not. bad) then
         do ip = 1, size(mh%pops)
            if (allocated(mh%pops(ip)%gsca) .neqv. allocated(mt%pops(ip)%gsca)) bad = .true.
         end do
      end if
      call fail_if(bad, 'which populations carry g', 0)
   end subroutine compare_routes


   subroutine report(name, a, b)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: a(:), b(:)
      real(wp) :: r, rmax
      integer  :: i, imax
      rmax = 0.0_wp;  imax = 1
      do i = 1, size(a)
         if (abs(b(i)) > 0.0_wp) then
            r = abs(a(i)/b(i) - 1.0_wp)
         else
            r = abs(a(i))
         end if
         if (r > rmax) then;  rmax = r;  imax = i;  end if
      end do
      write(*,'(a,a,a,es9.2,a,i0,a)') '   ', name, '  max rel ', rmax, &
           '  (at node ', imax, ')'
      call fail_if(rmax > TOL, name, 0)
   end subroutine report


   subroutine fail_if(cond, what, code)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: what
      integer,          intent(in) :: code
      if (.not. cond) return
      if (code /= 0) then
         write(*,'(a,a,a,i0)') '   *** ', what, ' = ', code
      else
         write(*,'(a,a)') '   *** over tolerance: ', what
      end if
      nbad = nbad + 1
   end subroutine fail_if

end program test_dustem_product
