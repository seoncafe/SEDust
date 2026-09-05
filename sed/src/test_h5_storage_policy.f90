program test_h5_storage_policy
   ! Does every HDF5 product in this tree obey the storage-precision policy?
   !
   !   ./test_h5_storage_policy.x
   !
   ! The policy, stated once in sed/src/sedust_h5.f90 above the writers and
   ! enforced here, has three classes decided by what a dataset IS:
   !
   !   quantities  Q_abs, Q_sca, Q_ext, Q_re, g, Z, F_tot, F_ref, the C_*
   !               cross sections, K_abs, albedo -- anything a calculation
   !               produced.  32-bit float in the file.  The physics carries
   !               its own uncertainty far above float32's 1.2e-7 relative
   !               resolution, so the other four bytes buy nothing.
   !
   !   axes        lambda, a_eff, theta, theta_i, theta_s, phi, band_lambda.
   !               64-bit float.  Their problem is DISTINCTNESS, not accuracy:
   !               the astrodust wavelength grid resolves each X-ray absorption
   !               edge with a pair of points 6.7e-7 apart in relative
   !               wavelength, only six to eleven representable float32 values
   !               apart, and a grid that stopped being strictly ascending has
   !               broken this code once already.
   !
   !   codes       regime.  One signed byte.
   !
   ! This test is what keeps an axis from silently becoming single later: the
   ! writers take `single` as an OPTIONAL argument defaulting to lossless, so a
   ! wrong tag is a one-word edit and nothing else would notice it.  It reads
   ! the FILE types only -- the memory type on both write and read is still
   ! double, which is why no reader had to change.
   use, intrinsic :: iso_fortran_env, only: int64
#ifdef SEDUST_HDF5
   use hdf5
#endif
   implicit none

   ! Every model this tree ships.  A model whose product is absent is skipped
   ! with a note rather than failed: not every tree carries every product.
   character(len=*), parameter :: MODELS(6) = &
      [character(len=9) :: 'astrodust', 'dl07', 'mrn', 'zubko', 'themis', 'g18d']

   ! The coordinate axes, matched on the dataset's own name because the same
   ! axis appears under several groups.
   character(len=*), parameter :: AXES(7) = &
      [character(len=11) :: 'lambda', 'a_eff', 'theta', 'theta_i', 'theta_s', &
                            'phi', 'band_lambda']

   integer :: nbad, nfile, nds

#ifndef SEDUST_HDF5
   write(*,'(a)') ' test_h5_storage_policy: built without HDF5; nothing to check'
   stop 0
#else
   integer :: e, i
   logical :: there
   character(len=256) :: path

   call h5open_f(e)
   if (e /= 0) then
      write(*,'(a)') ' test_h5_storage_policy: cannot start HDF5';  stop 1
   end if

   nbad = 0;  nfile = 0;  nds = 0
   do i = 1, size(MODELS)
      path = '../data/'//trim(MODELS(i))//'/sedust_'//trim(MODELS(i))//'.h5'
      inquire(file=trim(path), exist=there)
      if (.not. there) then
         write(*,'(a,a)') ' skip (absent): ', trim(path)
         cycle
      end if
      nfile = nfile + 1
      call check_file(trim(path))
   end do

   call h5close_f(e)

   write(*,'(a)') ''
   write(*,'(a,i0,a,i0,a)') ' checked ', nds, ' datasets in ', nfile, ' product(s)'
   if (nfile == 0) then
      write(*,'(a)') ' NO PRODUCTS FOUND';  stop 1
   end if
   if (nbad == 0) then
      write(*,'(a)') ' PASS'
   else
      write(*,'(i0,a)') nbad, ' DATASET(S) VIOLATE THE STORAGE POLICY'
      stop 1
   end if

contains

   subroutine check_file(path)
      character(len=*), intent(in) :: path
      integer(hid_t) :: fid
      integer :: e2, bad0
      bad0 = nbad
      call h5fopen_f(path, H5F_ACC_RDONLY_F, fid, e2)
      if (e2 /= 0) then
         write(*,'(a,a)') ' cannot open ', path
         nbad = nbad + 1;  return
      end if
      call walk(fid, '')
      call h5fclose_f(fid, e2)
      if (nbad == bad0) then
         write(*,'(a,a)') ' ok   ', path
      else
         write(*,'(a,a)') ' BAD  ', path
      end if
   end subroutine check_file


   recursive subroutine walk(loc, prefix)
      ! Every dataset under loc, depth first.
      integer(hid_t),   intent(in) :: loc
      character(len=*), intent(in) :: prefix
      integer(hid_t) :: gid
      integer :: e2, n, k, otype, ln
      character(len=128) :: nm
      call h5gn_members_f(loc, '.', n, e2)
      if (e2 /= 0) return
      do k = 0, n-1
         call h5gget_obj_info_idx_f(loc, '.', k, nm, otype, e2)
         if (e2 /= 0) cycle
         ln = len_trim(nm)
         if (otype == H5G_GROUP_F) then
            call h5gopen_f(loc, nm(1:ln), gid, e2)
            if (e2 == 0) then
               call walk(gid, prefix//nm(1:ln)//'/')
               call h5gclose_f(gid, e2)
            end if
         else if (otype == H5G_DATASET_F) then
            call check_dataset(loc, nm(1:ln), prefix//nm(1:ln))
         end if
      end do
   end subroutine walk


   subroutine check_dataset(loc, name, full)
      integer(hid_t),   intent(in) :: loc
      character(len=*), intent(in) :: name, full
      integer(hid_t)  :: did, tid
      integer(size_t) :: nbyte
      integer :: e2, tclass, want_class, want_byte
      character(len=16) :: what
      call h5dopen_f(loc, name, did, e2)
      if (e2 /= 0) then
         write(*,'(a,a)') '   cannot open dataset ', full
         nbad = nbad + 1;  return
      end if
      call h5dget_type_f(did, tid, e2)
      call h5tget_class_f(tid, tclass, e2)
      call h5tget_size_f(tid, nbyte, e2)
      call h5tclose_f(tid, e2)
      call h5dclose_f(did, e2)
      nds = nds + 1

      if (name == 'regime') then
         what = 'code';       want_class = H5T_INTEGER_F;  want_byte = 1
      else if (is_axis(name)) then
         what = 'axis';       want_class = H5T_FLOAT_F;    want_byte = 8
      else
         what = 'quantity';   want_class = H5T_FLOAT_F;    want_byte = 4
      end if

      if (tclass /= want_class .or. int(nbyte) /= want_byte) then
         write(*,'(a,a,a,a,a,i0,a,a,a,i0,a)') &
            '   ', full, ' (', trim(what), '): file type is ', &
            8*int(nbyte), '-bit ', trim(class_name(tclass)), &
            ', policy says ', 8*want_byte, '-bit '//trim(class_name(want_class))
         nbad = nbad + 1
      end if
   end subroutine check_dataset


   logical function is_axis(name)
      character(len=*), intent(in) :: name
      integer :: k
      is_axis = .false.
      do k = 1, size(AXES)
         if (name == trim(AXES(k))) then
            is_axis = .true.;  return
         end if
      end do
   end function is_axis


   function class_name(c) result(s)
      integer, intent(in) :: c
      character(len=16) :: s
      if (c == H5T_FLOAT_F) then
         s = 'float'
      else if (c == H5T_INTEGER_F) then
         s = 'integer'
      else
         s = 'other'
      end if
   end function class_name
#endif

end program test_h5_storage_policy
