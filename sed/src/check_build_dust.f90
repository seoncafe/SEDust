program check_build_dust
   ! Does build_dust, reading the HDF5 products, give the same model as the
   ! builders reading the text products?
   !
   !   ./check_build_dust.x [astrodust | dl07 | zubko | all]
   !
   ! For each model and each wavelength set, this builds the model twice --
   ! once through build_dust on data/<model>/sedust_<model>.h5 and once through the
   ! model's own builder on the text files -- and compares the grid, the size
   ! integral (C_ext, C_abs, C_sca, <cos>), the dust mass per H, and the
   ! extinction curve dust_extinction serves.
   !
   ! The two are not required to agree bit for bit and are not expected to: the
   ! text tables carry seven significant digits and the HDF5 carries the full
   ! double the optics were computed in, so the HDF5 route is the more exact of
   ! the two and the difference between them is the text precision.  What is
   ! checked is that nothing is off by more than that.
   use constants,         only: wp
   use sed_astrodust_mod, only: dust_model_t, build_dust, build_astrodust, build_dl07, &
                                build_zubko, size_integrated_extinction, &
                                dust_mass_per_H, dust_extinction, &
                                dl07_euv_lambda_floor, sed_verbose
   use euv_astrodust_tmatrix, only: use_tmatrix_euv_band_optics
   use zubko_io,          only: read_zubko_optics
   use sedust_product_mod, only: read_sedust_qtable
   implicit none

   character(len=*), parameter :: DDIR   = '../data'
   ! The TEXT side of the comparison.  Everything else in the tree now reads the
   ! HDF5 product, which is exactly what this program checks the text against,
   ! so these two are among the last readers of the text tables.  They are not
   ! shipped: run sed/calc_qtable.x (and the T-matrix sweep for the astrodust
   ! pair) before this test.
   character(len=*), parameter :: F_QT   = '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400.dat'
   character(len=*), parameter :: F_QT_E = '../data/astrodust/q_astrodust_P0.20_Fe0.00_1.400_euv.dat'
   character(len=*), parameter :: F_SD   = '../data/release/size_distribution.dat'
   character(len=*), parameter :: F_CFG  = '../data/zubko/ZDA_BARE_GR_S_Config.dat'
   character(len=*), parameter :: D_ZUB  = '../data/zubko/'
   character(len=*), parameter :: F_ZU_H5 = '../data/zubko/sedust_zubko.h5'
   character(len=*), parameter :: K_AD   = '../data/astrodust/kext_astrodust_MW_euv.dat'
   character(len=*), parameter :: K_DL   = '../data/dl07/kext_dl07_MW_euv.dat'
   ! The mie_d03 curve, because the text side of the zubko comparison is the
   ! mie_d03 optics -- the text q_zubko_*.dat ARE that set.  Naming the default
   ! curve here would compare the size integral of one set of optics against
   ! the stored integral of the other.
   character(len=*), parameter :: K_ZU   = &
        '../data/zubko/kext_zubko_BARE_GR_S_mie_d03_euv.dat'
   integer,  parameter :: NT_IN = 100
   real(wp), parameter :: T_LO = 1.0_wp, T_HI = 3000.0_wp
   ! Seven written digits in the text tables: a value can be off by 5e-7 of
   ! itself, and a ratio of two such by twice that.  Stand a factor of 10 above
   ! it, so that this reports a real disagreement and not the last digit.
   real(wp), parameter :: TOL = 1.0e-5_wp

   character(len=32) :: which
   integer :: narg, nbad

   call use_tmatrix_euv_band_optics()
   sed_verbose = .false.
   narg = command_argument_count()
   which = 'all'
   if (narg >= 1) call get_command_argument(1, which)

   nbad = 0
   select case (trim(which))
   case ('astrodust');  call check_astrodust()
   case ('dl07');       call check_dl07()
   case ('zubko');      call check_zubko()
   case ('all');        call check_astrodust();  call check_dl07();  call check_zubko()
   case default
      write(*,'(a)') ' usage: ./check_build_dust.x [astrodust | dl07 | zubko | all]'
      stop 1
   end select

   write(*,'(a)') ''
   if (nbad == 0) then
      write(*,'(a)') ' ALL CHECKS PASSED'
   else
      write(*,'(i0,a)') nbad, ' CHECK(S) FAILED'
      stop 1
   end if

contains

   subroutine check_astrodust()
      type(dust_model_t) :: mh, mt
      integer :: st
      call build_dust(mh, 'astrodust', DDIR, NT_IN, T_LO, T_HI, .false., status=st)
      call fail_if(st /= 0, 'astrodust narrow: build_dust status', st)
      call build_astrodust(mt, F_QT, F_SD, NT_IN, T_LO, T_HI, status=st, kext_path=K_AD)
      call fail_if(st /= 0, 'astrodust narrow: build_astrodust status', st)
      call compare('astrodust, non-EUV', mh, mt)

      call build_dust(mh, 'astrodust', DDIR, NT_IN, T_LO, T_HI, .true., status=st)
      call fail_if(st /= 0, 'astrodust wide: build_dust status', st)
      call build_astrodust(mt, F_QT_E, F_SD, NT_IN, T_LO, T_HI, status=st, kext_path=K_AD)
      call fail_if(st /= 0, 'astrodust wide: build_astrodust status', st)
      call compare('astrodust, EUV', mh, mt)
   end subroutine check_astrodust

   subroutine check_dl07()
      type(dust_model_t) :: mh, mt
      integer :: st
      call build_dust(mh, 'dl07', DDIR, NT_IN, T_LO, T_HI, .false., status=st)
      call fail_if(st /= 0, 'dl07 narrow: build_dust status', st)
      call build_dl07(mt, F_QT, F_SD, 7, 1.0_wp, NT_IN, T_LO, T_HI, status=st, kext_path=K_DL)
      call fail_if(st /= 0, 'dl07 narrow: build_dl07 status', st)
      call compare('dl07, non-EUV', mh, mt)

      call build_dust(mh, 'dl07', DDIR, NT_IN, T_LO, T_HI, .true., status=st)
      call fail_if(st /= 0, 'dl07 wide: build_dust status', st)
      call build_dl07(mt, F_QT_E, F_SD, 7, 1.0_wp, NT_IN, T_LO, T_HI, status=st, &
                      lam_min=dl07_euv_lambda_floor(), kext_path=K_DL)
      call fail_if(st /= 0, 'dl07 wide: build_dl07 status', st)
      call compare('dl07, EUV', mh, mt)
   end subroutine check_dl07

   subroutine check_zubko()
      ! zubko stores TWO optics sets, so this compares like with like: the text
      ! q_zubko_*.dat products ARE the mie_d03 recomputation, so that is the
      ! set the HDF5 side is asked for.  The DEFAULT set is the benchmark's own
      ! tables, which have no text counterpart to compare against here -- they
      ! ARE text, and check_zubko_zda below reads them back directly.
      type(dust_model_t) :: mh, mt
      integer :: st
      call build_dust(mh, 'zubko', DDIR, NT_IN, T_LO, T_HI, .false., status=st, &
                      zubko_optics='mie_d03')
      call fail_if(st /= 0, 'zubko narrow: build_dust status', st)
      call build_zubko(mt, F_CFG, D_ZUB, NT_IN, T_LO, T_HI, status=st, kext_path=K_ZU, &
                       include_euv=.false.)
      call fail_if(st /= 0, 'zubko narrow: build_zubko status', st)
      call compare('zubko mie_d03, non-EUV', mh, mt)

      call build_dust(mh, 'zubko', DDIR, NT_IN, T_LO, T_HI, .true., status=st, &
                      zubko_optics='mie_d03')
      call fail_if(st /= 0, 'zubko wide: build_dust status', st)
      call build_zubko(mt, F_CFG, D_ZUB, NT_IN, T_LO, T_HI, status=st, kext_path=K_ZU, &
                       include_euv=.true.)
      call fail_if(st /= 0, 'zubko wide: build_zubko status', st)
      call compare('zubko mie_d03, EUV', mh, mt)

      call check_zubko_zda()
   end subroutine check_zubko


   subroutine check_zubko_zda()
      ! Does /qtable/{sil,gra,pah} hold the distributed tables as they stand?
      ! This is the set the model is built on by default, and the seven codes
      ! the SHG benchmark compares against read those same files, so anything
      ! but an exact round-trip would make that comparison measure an optics
      ! difference instead of the solver.  Read both sides and diff them.
      character(len=*), parameter :: COMP(3) = [character(len=3) :: 'sil', 'gra', 'pah']
      character(len=*), parameter :: FILE(3) = [character(len=24) :: &
           'suvSil_121_1201.dat', 'Gra_121_1201.dat', 'PAH_28_1201_neu.dat']
      real(wp), allocatable :: a_f(:), l_f(:), qa_f(:,:), qs_f(:,:), gg_f(:,:)
      real(wp), allocatable :: a_h(:), qe_h(:,:), qa_h(:,:), qs_h(:,:), gg_h(:,:)
      real(wp) :: rho_f, rho_h, da, ds, dg
      integer  :: ic, nsize, nwave
      logical  :: ok

      write(*,'(a)') ''
      write(*,'(a)') ' === zubko default set: /qtable vs the distributed files'
      do ic = 1, 3
         call read_zubko_optics(D_ZUB//trim(FILE(ic)), nsize, nwave, a_f, l_f, &
                                qa_f, qs_f, rho_f, ok=ok, gpar=gg_f)
         call fail_if(.not. ok, 'zubko zda: cannot read '//trim(FILE(ic)), 0)
         call read_sedust_qtable(F_ZU_H5, trim(COMP(ic)), .true., a_h, qe_h, qa_h, &
                                 qs_h, gg_h, rho_h, ok)
         call fail_if(.not. ok, 'zubko zda: cannot read /qtable/'//trim(COMP(ic)), 0)
         if (size(a_h) /= nsize .or. size(qa_h,1) /= nwave) then
            write(*,'(a,a)') '   *** shape mismatch for ', trim(COMP(ic))
            nbad = nbad + 1;  cycle
         end if
         da = maxval(abs(qa_h - qa_f));  ds = maxval(abs(qs_h - qs_f))
         dg = maxval(abs(gg_h - gg_f))
         write(*,'(a,a4,a,es9.2,a,es9.2,a,es9.2)') '   ', trim(COMP(ic)), &
              '  max|dQ_abs| ', da, '  max|dQ_sca| ', ds, '  max|dg| ', dg
         ! Stored as written: the product holds these arrays, it does not
         ! recompute them, so anything but zero is a transcription bug.
         call fail_if(da /= 0.0_wp .or. ds /= 0.0_wp .or. dg /= 0.0_wp, &
                      'zubko zda: /qtable/'//trim(COMP(ic))//' is not the file', 0)
         deallocate(a_f, l_f, qa_f, qs_f, gg_f, a_h, qe_h, qa_h, qs_h, gg_h)
      end do
   end subroutine check_zubko_zda


   subroutine compare(label, mh, mt)
      character(len=*),   intent(in) :: label
      type(dust_model_t), intent(in) :: mh, mt
      real(wp), allocatable :: eh(:), ah(:), sh(:), gh(:)
      real(wp), allocatable :: et(:), at(:), st_(:), gt(:)
      real(wp), allocatable :: xe(:), xa(:), xs(:), xg(:)
      real(wp), allocatable :: ye(:), ya(:), ys(:), yg(:)
      real(wp) :: dlam, mdh, mdt
      integer  :: n, s1, s2

      write(*,'(a)') ''
      write(*,'(a,a)') ' === ', label
      if (mh%NLAM /= mt%NLAM) then
         write(*,'(a,i0,a,i0)') '   NLAM differs: HDF5 ', mh%NLAM, ', text ', mt%NLAM
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

      mdh = dust_mass_per_H(mh);  mdt = dust_mass_per_H(mt)
      write(*,'(a,es14.7,a,es14.7,a,es9.2)') '   M_dust/H  HDF5 ', mdh, '  text ', mdt, &
         '   rel ', abs(mdh/mdt - 1.0_wp)
      call fail_if(abs(mdh/mdt - 1.0_wp) > TOL, 'M_dust/H', 0)

      ! The curve dust_extinction serves: from /kext on one side and from the
      ! text table on the other, interpolated onto the same model grid.
      allocate(xe(n), xa(n), xs(n), xg(n), ye(n), ya(n), ys(n), yg(n))
      call dust_extinction(mh, xe, xa, xs, gbar=xg, status=s1)
      call dust_extinction(mt, ye, ya, ys, gbar=yg, status=s2)
      if (s1 /= 0 .or. s2 /= 0) then
         write(*,'(a,i0,a,i0)') '   *** dust_extinction status: HDF5 ', s1, ', text ', s2
         nbad = nbad + 1
      else
         call report('kext C_ext', xe, ye)
         call report('kext C_abs', xa, ya)
         call report('kext <cos>', xg, yg)
      end if
      deallocate(eh, ah, sh, gh, et, at, st_, gt, xe, xa, xs, xg, ye, ya, ys, yg)
   end subroutine compare


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
      write(*,'(a,a,a,es9.2,a,i0,a)') '   ', name, '  max rel ', rmax, '  (at node ', imax, ')'
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

end program check_build_dust
