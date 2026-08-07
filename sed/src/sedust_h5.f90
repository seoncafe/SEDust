module sedust_h5
   ! The HDF5 layer under SEDust's optics products: the size-integrated
   ! extinction curves, the (lambda, a_eff) cross-section tables, and -- in the
   ! polarized branch -- the orientation-resolved and scattering-matrix
   ! products.  One file per dust model,
   !
   !   data/astrodust/sedust_astrodust.h5   data/dl07/sedust_dl07.h5   data/zubko/sedust_zubko.h5
   !
   ! because the models do not share a wavelength grid, a radius grid or a
   ! component list, and because a host ships only the models it runs.
   !
   ! WHY THE FILES CARRY THE IONIZING BAND AND THE HOST ASKS FOR IT.  Each file
   ! holds ONE wavelength axis, the widest the model has, and records the index
   ! at which the non-ionizing part begins:
   !
   !   /grid/lambda      (n_lam)  [um], strictly increasing
   !   /grid  attribute  i_lyman  = first index with lambda >= 0.0912 um
   !
   ! A reader asked for include_euv = .false. returns lambda(i_lyman:) and the
   ! same slice of every wavelength-indexed array.  The narrow product is then
   ! a VIEW of the wide one rather than a second file that has to be kept in
   ! step with it: they cannot disagree, and there is nothing to regenerate
   ! twice.  This replaces the paired *_euv.dat text products.
   !
   ! ARRAY SHAPES.  Fortran declares cross sections (n_lam, n_a) -- the order
   ! the solver already uses -- so on disk the dimensions are [n_lam, n_a] in
   ! HDF5's Fortran convention, which h5py reads back as (n_a, n_lam).  Each
   ! grain radius is then one contiguous spectrum in Python, and the EUV cut is
   ! a hyperslab in the fastest-varying Fortran dimension.
   !
   ! PROVENANCE IS PART OF THE FORMAT.  Every group carries attributes naming
   ! what produced it and from which inputs.  That is not decoration: the ZDA
   ! optics tables this tree reads state a dielectric function in their header
   ! that they were demonstrably NOT computed from (see zubko_io.f90), and a
   ! stale label of that kind costs real time to unpick.  Writers here record
   ! the file paths they actually opened.
   !
   ! Text products are still written alongside; HDF5 is the primary form and
   ! the text is for reading with an eye and diffing.
   use, intrinsic :: iso_fortran_env, only: real64, int64
#ifdef SEDUST_HDF5
   use hdf5
#endif
   implicit none
   private

   ! Handle kind.  With HDF5 compiled in it IS hid_t, so every call below type
   ! checks against the library; without it the module still compiles and every
   ! entry point reports "not available" through its ok argument, which is what
   ! the callers already do when a file is missing.  Nothing outside needs an
   ! #ifdef of its own.
#ifdef SEDUST_HDF5
   integer, parameter, public :: h5id_k = hid_t
   logical, parameter, public :: sedust_has_hdf5 = .true.
#else
   integer, parameter, public :: h5id_k = int64
   logical, parameter, public :: sedust_has_hdf5 = .false.
#endif

   public :: h5_begin, h5_end
   public :: h5_create_file, h5_open_file, h5_open_rw, h5_close_file
   public :: h5_group, h5_group_open, h5_group_close, h5_unlink
   public :: h5_write_1d, h5_write_2d, h5_write_3d, h5_write_6d
   public :: h5_read_1d, h5_read_2d, h5_read_3d, h5_read_6d
   public :: h5_dims
   public :: h5_put_attr_d, h5_put_attr_i, h5_put_attr_s
   public :: h5_get_attr_d, h5_get_attr_i, h5_get_attr_s
   public :: h5_has

#ifdef SEDUST_HDF5
   ! HDF5 must be opened once per process and closed once.  Reference-counted
   ! so that independent callers (a writer, then a reader) do not close the
   ! library out from under each other.
   integer, save :: n_open = 0
#endif

   ! Chunking and compression.  A chunk spans the wavelength axis in blocks so
   ! that an EUV cut reads whole chunks, and the whole radius axis so that one
   ! grain's spectrum never straddles chunks.
   integer, parameter :: LAM_CHUNK = 256
   integer, parameter :: GZIP_LEVEL = 4

contains

#ifdef SEDUST_HDF5

   subroutine h5_begin(ok)
      logical, intent(out) :: ok
      integer :: e
      ok = .true.
      if (n_open == 0) then
         call h5open_f(e)
         ok = (e == 0)
         if (.not. ok) return
      end if
      n_open = n_open + 1
   end subroutine h5_begin

   subroutine h5_end()
      integer :: e
      if (n_open <= 0) return
      n_open = n_open - 1
      if (n_open == 0) call h5close_f(e)
   end subroutine h5_end


   subroutine h5_create_file(path, fid, ok)
      character(len=*), intent(in)  :: path
      integer(h5id_k),   intent(out) :: fid
      logical,          intent(out) :: ok
      integer :: e
      call h5fcreate_f(trim(path), H5F_ACC_TRUNC_F, fid, e)
      ok = (e == 0)
   end subroutine h5_create_file

   subroutine h5_open_file(path, fid, ok)
      character(len=*), intent(in)  :: path
      integer(h5id_k),   intent(out) :: fid
      logical,          intent(out) :: ok
      integer :: e
      call h5fopen_f(trim(path), H5F_ACC_RDONLY_F, fid, e)
      ok = (e == 0)
   end subroutine h5_open_file

   subroutine h5_open_rw(path, fid, ok)
      ! Open an existing file for appending -- the writers add one component
      ! group at a time rather than holding the whole model in memory.
      character(len=*), intent(in)  :: path
      integer(h5id_k),  intent(out) :: fid
      logical,          intent(out) :: ok
      integer :: e
      call h5fopen_f(trim(path), H5F_ACC_RDWR_F, fid, e)
      ok = (e == 0)
   end subroutine h5_open_rw

   subroutine h5_close_file(fid)
      integer(h5id_k), intent(in) :: fid
      integer :: e
      call h5fclose_f(fid, e)
   end subroutine h5_close_file


   subroutine h5_group(loc, name, gid, ok)
      ! Create a group (intermediate groups are created too).
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      integer(h5id_k),   intent(out) :: gid
      logical,          intent(out) :: ok
      integer(h5id_k) :: lcpl
      integer :: e
      call h5pcreate_f(H5P_LINK_CREATE_F, lcpl, e)
      call h5pset_create_inter_group_f(lcpl, 1, e)
      call h5gcreate_f(loc, trim(name), gid, e, lcpl_id=lcpl)
      call h5pclose_f(lcpl, e)
      ok = (e == 0)
   end subroutine h5_group

   subroutine h5_group_open(loc, name, gid, ok)
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      integer(h5id_k),   intent(out) :: gid
      logical,          intent(out) :: ok
      integer :: e
      call h5gopen_f(loc, trim(name), gid, e)
      ok = (e == 0)
   end subroutine h5_group_open

   subroutine h5_group_close(gid)
      integer(h5id_k), intent(in) :: gid
      integer :: e
      call h5gclose_f(gid, e)
   end subroutine h5_group_close

   subroutine h5_unlink(loc, name)
      ! Drop a group or dataset so that a writer run twice replaces what it
      ! wrote rather than failing on the second pass.  HDF5 does not reclaim
      ! the space, so a file rewritten many times is worth regenerating; the
      ! writers here rewrite one group, not the file.
      integer(h5id_k),  intent(in) :: loc
      character(len=*), intent(in) :: name
      integer :: e
      call h5ldelete_f(loc, trim(name), e)
   end subroutine h5_unlink


   logical function h5_has(loc, name)
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name
      integer :: e
      call h5lexists_f(loc, trim(name), h5_has, e)
      if (e /= 0) h5_has = .false.
   end function h5_has


   ! ---- writers -------------------------------------------------------
   subroutine h5_write_1d(loc, name, a, units, long_name)
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name
      real(real64),     intent(in) :: a(:)
      character(len=*), intent(in), optional :: units, long_name
      integer(hsize_t) :: d(1)
      d = [int(size(a), hsize_t)]
      call write_real(loc, name, d, reshape(a, [size(a)]), units, long_name)
   end subroutine h5_write_1d

   subroutine h5_write_2d(loc, name, a, units, long_name)
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name
      real(real64),     intent(in) :: a(:,:)
      character(len=*), intent(in), optional :: units, long_name
      integer(hsize_t) :: d(2)
      d = [int(size(a,1), hsize_t), int(size(a,2), hsize_t)]
      call write_real(loc, name, d, reshape(a, [size(a)]), units, long_name)
   end subroutine h5_write_2d

   subroutine h5_write_3d(loc, name, a, units, long_name)
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name
      real(real64),     intent(in) :: a(:,:,:)
      character(len=*), intent(in), optional :: units, long_name
      integer(hsize_t) :: d(3)
      d = [int(size(a,1), hsize_t), int(size(a,2), hsize_t), int(size(a,3), hsize_t)]
      call write_real(loc, name, d, reshape(a, [size(a)]), units, long_name)
   end subroutine h5_write_3d

   subroutine h5_write_6d(loc, name, a, units, long_name)
      ! The aligned phase matrix, (n_ti, n_ts, n_phi, 4, 4, n_band).  Its chunk
      ! is ONE band: a host reads a band at a time, and the whole array as a
      ! single chunk would be tens of megabytes through the chunk cache for
      ! every such read.
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name
      real(real64),     intent(in) :: a(:,:,:,:,:,:)
      character(len=*), intent(in), optional :: units, long_name
      integer(hsize_t) :: d(6), ch(6)
      integer :: k
      do k = 1, 6
         d(k) = int(size(a,k), hsize_t)
      end do
      ch = d;  ch(6) = 1_hsize_t
      call write_real(loc, name, d, reshape(a, [size(a)]), units, long_name, chunk=ch)
   end subroutine h5_write_6d

   subroutine write_real(loc, name, d, flat, units, long_name, chunk)
      ! One writer for every rank: chunked and deflated, with the units and a
      ! one-line description attached to the dataset rather than to a README
      ! that can drift away from it.
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name
      integer(hsize_t), intent(in) :: d(:)
      real(real64),     intent(in) :: flat(:)
      character(len=*), intent(in), optional :: units, long_name
      ! Chunk shape, when the default -- blocks along the first (wavelength)
      ! axis and the whole of the others -- is the wrong one for this array.
      integer(hsize_t), intent(in), optional :: chunk(:)
      integer(h5id_k)   :: sid, did, dcpl
      integer(hsize_t) :: ch(size(d))
      integer :: e, r
      r = size(d)
      call h5screate_simple_f(r, d, sid, e)
      ch = d
      if (ch(1) > int(LAM_CHUNK, hsize_t)) ch(1) = int(LAM_CHUNK, hsize_t)
      if (present(chunk)) ch = chunk
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, e)
      call h5pset_chunk_f(dcpl, r, ch, e)
      call h5pset_shuffle_f(dcpl, e)
      call h5pset_deflate_f(dcpl, GZIP_LEVEL, e)
      call h5dcreate_f(loc, trim(name), H5T_NATIVE_DOUBLE, sid, did, e, dcpl_id=dcpl)
      call h5dwrite_f(did, H5T_NATIVE_DOUBLE, flat, d, e)
      if (present(units))     call h5_put_attr_s(did, 'units', units)
      if (present(long_name)) call h5_put_attr_s(did, 'long_name', long_name)
      call h5dclose_f(did, e)
      call h5pclose_f(dcpl, e)
      call h5sclose_f(sid, e)
   end subroutine write_real


   ! ---- readers -------------------------------------------------------
   subroutine h5_dims(loc, name, d, ok)
      ! Shape of a dataset, so a caller can allocate before reading.
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      integer(hsize_t), allocatable, intent(out) :: d(:)
      logical,          intent(out) :: ok
      integer(h5id_k)   :: did, sid
      integer(hsize_t) :: dm(7), mx(7)
      integer :: e, r
      ok = .false.
      call h5dopen_f(loc, trim(name), did, e);  if (e /= 0) return
      call h5dget_space_f(did, sid, e)
      call h5sget_simple_extent_ndims_f(sid, r, e)
      call h5sget_simple_extent_dims_f(sid, dm, mx, e)
      allocate(d(r));  d = dm(1:r)
      call h5sclose_f(sid, e);  call h5dclose_f(did, e)
      ok = .true.
   end subroutine h5_dims

   subroutine h5_read_1d(loc, name, a, ok, i0)
      ! Read a 1-D dataset, optionally from index i0 to the end -- the EUV cut.
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      real(real64), allocatable, intent(out) :: a(:)
      logical,          intent(out) :: ok
      integer, optional, intent(in) :: i0
      integer(hsize_t), allocatable :: d(:)
      integer(hsize_t) :: off(1), cnt(1)
      integer(h5id_k)   :: did, sid, msid
      integer :: e, lo
      ok = .false.
      call h5_dims(loc, name, d, ok);  if (.not. ok) return
      lo = 1;  if (present(i0)) lo = max(1, i0)
      cnt = [d(1) - int(lo-1, hsize_t)]
      allocate(a(cnt(1)))
      call h5dopen_f(loc, trim(name), did, e)
      call h5dget_space_f(did, sid, e)
      off = [int(lo-1, hsize_t)]
      call h5sselect_hyperslab_f(sid, H5S_SELECT_SET_F, off, cnt, e)
      call h5screate_simple_f(1, cnt, msid, e)
      call h5dread_f(did, H5T_NATIVE_DOUBLE, a, cnt, e, mem_space_id=msid, file_space_id=sid)
      call h5sclose_f(msid, e);  call h5sclose_f(sid, e);  call h5dclose_f(did, e)
      ok = (e == 0)
   end subroutine h5_read_1d

   subroutine h5_read_2d(loc, name, a, ok, i0)
      ! (n_lam, n_a); i0 cuts the FIRST dimension, which is the wavelength.
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      real(real64), allocatable, intent(out) :: a(:,:)
      logical,          intent(out) :: ok
      integer, optional, intent(in) :: i0
      integer(hsize_t), allocatable :: d(:)
      integer(hsize_t) :: off(2), cnt(2)
      integer(h5id_k)   :: did, sid, msid
      integer :: e, lo
      ok = .false.
      call h5_dims(loc, name, d, ok);  if (.not. ok) return
      lo = 1;  if (present(i0)) lo = max(1, i0)
      cnt = [d(1) - int(lo-1, hsize_t), d(2)]
      allocate(a(cnt(1), cnt(2)))
      call h5dopen_f(loc, trim(name), did, e)
      call h5dget_space_f(did, sid, e)
      off = [int(lo-1, hsize_t), 0_hsize_t]
      call h5sselect_hyperslab_f(sid, H5S_SELECT_SET_F, off, cnt, e)
      call h5screate_simple_f(2, cnt, msid, e)
      call h5dread_f(did, H5T_NATIVE_DOUBLE, a, cnt, e, mem_space_id=msid, file_space_id=sid)
      call h5sclose_f(msid, e);  call h5sclose_f(sid, e);  call h5dclose_f(did, e)
      ok = (e == 0)
   end subroutine h5_read_2d

   subroutine h5_read_3d(loc, name, a, ok, i0)
      ! (n_lam, n_a, n_k); i0 cuts the wavelength dimension as above.
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      real(real64), allocatable, intent(out) :: a(:,:,:)
      logical,          intent(out) :: ok
      integer, optional, intent(in) :: i0
      integer(hsize_t), allocatable :: d(:)
      integer(hsize_t) :: off(3), cnt(3)
      integer(h5id_k)   :: did, sid, msid
      integer :: e, lo
      ok = .false.
      call h5_dims(loc, name, d, ok);  if (.not. ok) return
      lo = 1;  if (present(i0)) lo = max(1, i0)
      cnt = [d(1) - int(lo-1, hsize_t), d(2), d(3)]
      allocate(a(cnt(1), cnt(2), cnt(3)))
      call h5dopen_f(loc, trim(name), did, e)
      call h5dget_space_f(did, sid, e)
      off = [int(lo-1, hsize_t), 0_hsize_t, 0_hsize_t]
      call h5sselect_hyperslab_f(sid, H5S_SELECT_SET_F, off, cnt, e)
      call h5screate_simple_f(3, cnt, msid, e)
      call h5dread_f(did, H5T_NATIVE_DOUBLE, a, cnt, e, mem_space_id=msid, file_space_id=sid)
      call h5sclose_f(msid, e);  call h5sclose_f(sid, e);  call h5dclose_f(did, e)
      ok = (e == 0)
   end subroutine h5_read_3d

   subroutine h5_read_6d(loc, name, a, ok)
      ! The aligned phase matrix, (n_ti, n_ts, n_phi, 4, 4, n_band).  No i0:
      ! this array has no wavelength axis to cut -- its last dimension is the
      ! handful of bands the aligned table was computed in.
      integer(h5id_k),  intent(in)  :: loc
      character(len=*), intent(in)  :: name
      real(real64), allocatable, intent(out) :: a(:,:,:,:,:,:)
      logical,          intent(out) :: ok
      integer(hsize_t), allocatable :: d(:)
      integer(hsize_t) :: cnt(6)
      integer(h5id_k)  :: did
      integer :: e, k
      ok = .false.
      call h5_dims(loc, name, d, ok);  if (.not. ok) return
      if (size(d) /= 6) then;  ok = .false.;  return;  end if
      do k = 1, 6
         cnt(k) = d(k)
      end do
      allocate(a(cnt(1), cnt(2), cnt(3), cnt(4), cnt(5), cnt(6)))
      call h5dopen_f(loc, trim(name), did, e)
      call h5dread_f(did, H5T_NATIVE_DOUBLE, a, cnt, e)
      call h5dclose_f(did, e)
      ok = (e == 0)
   end subroutine h5_read_6d


   ! ---- attributes ----------------------------------------------------
   subroutine h5_put_attr_d(loc, name, v)
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name
      real(real64),     intent(in) :: v
      integer(h5id_k) :: sid, aid
      integer(hsize_t) :: d(1)
      integer :: e
      d = [1_hsize_t]
      call h5screate_f(H5S_SCALAR_F, sid, e)
      call h5acreate_f(loc, trim(name), H5T_NATIVE_DOUBLE, sid, aid, e)
      call h5awrite_f(aid, H5T_NATIVE_DOUBLE, v, d, e)
      call h5aclose_f(aid, e);  call h5sclose_f(sid, e)
   end subroutine h5_put_attr_d

   subroutine h5_put_attr_i(loc, name, v)
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name
      integer,          intent(in) :: v
      integer(h5id_k) :: sid, aid
      integer(hsize_t) :: d(1)
      integer :: e
      d = [1_hsize_t]
      call h5screate_f(H5S_SCALAR_F, sid, e)
      call h5acreate_f(loc, trim(name), H5T_NATIVE_INTEGER, sid, aid, e)
      call h5awrite_f(aid, H5T_NATIVE_INTEGER, v, d, e)
      call h5aclose_f(aid, e);  call h5sclose_f(sid, e)
   end subroutine h5_put_attr_i

   subroutine h5_put_attr_s(loc, name, v)
      integer(h5id_k),   intent(in) :: loc
      character(len=*), intent(in) :: name, v
      integer(h5id_k) :: sid, aid, tid
      integer(hsize_t) :: d(1)
      integer :: e
      d = [1_hsize_t]
      call h5tcopy_f(H5T_NATIVE_CHARACTER, tid, e)
      call h5tset_size_f(tid, int(max(len_trim(v),1), size_t), e)
      call h5screate_f(H5S_SCALAR_F, sid, e)
      call h5acreate_f(loc, trim(name), tid, sid, aid, e)
      call h5awrite_f(aid, tid, trim(v), d, e)
      call h5aclose_f(aid, e);  call h5sclose_f(sid, e);  call h5tclose_f(tid, e)
   end subroutine h5_put_attr_s

   subroutine h5_get_attr_d(loc, name, v, ok)
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      real(real64),     intent(out) :: v
      logical,          intent(out) :: ok
      integer(h5id_k) :: aid
      integer(hsize_t) :: d(1)
      integer :: e
      d = [1_hsize_t];  ok = .false.;  v = 0.0_real64
      call h5aopen_f(loc, trim(name), aid, e);  if (e /= 0) return
      call h5aread_f(aid, H5T_NATIVE_DOUBLE, v, d, e)
      call h5aclose_f(aid, e);  ok = (e == 0)
   end subroutine h5_get_attr_d

   subroutine h5_get_attr_i(loc, name, v, ok)
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      integer,          intent(out) :: v
      logical,          intent(out) :: ok
      integer(h5id_k) :: aid
      integer(hsize_t) :: d(1)
      integer :: e
      d = [1_hsize_t];  ok = .false.;  v = 0
      call h5aopen_f(loc, trim(name), aid, e);  if (e /= 0) return
      call h5aread_f(aid, H5T_NATIVE_INTEGER, v, d, e)
      call h5aclose_f(aid, e);  ok = (e == 0)
   end subroutine h5_get_attr_i

   subroutine h5_get_attr_s(loc, name, v, ok)
      integer(h5id_k),   intent(in)  :: loc
      character(len=*), intent(in)  :: name
      character(len=*), intent(out) :: v
      logical,          intent(out) :: ok
      integer(h5id_k) :: aid, tid
      integer(hsize_t) :: d(1)
      integer :: e
      d = [1_hsize_t];  ok = .false.;  v = ''
      call h5aopen_f(loc, trim(name), aid, e);  if (e /= 0) return
      call h5aget_type_f(aid, tid, e)
      call h5aread_f(aid, tid, v, d, e)
      call h5tclose_f(tid, e);  call h5aclose_f(aid, e)
      ok = (e == 0)
   end subroutine h5_get_attr_s

#else
   ! ---- built without HDF5: every entry point declines ------------------
   subroutine h5_begin(ok);  logical, intent(out) :: ok;  ok = .false.;  end subroutine
   subroutine h5_end();      end subroutine
   subroutine h5_create_file(path, fid, ok)
      character(len=*), intent(in) :: path;  integer(h5id_k), intent(out) :: fid
      logical, intent(out) :: ok;  fid = 0;  ok = .false.
   end subroutine
   subroutine h5_open_file(path, fid, ok)
      character(len=*), intent(in) :: path;  integer(h5id_k), intent(out) :: fid
      logical, intent(out) :: ok;  fid = 0;  ok = .false.
   end subroutine
   subroutine h5_open_rw(path, fid, ok)
      character(len=*), intent(in) :: path;  integer(h5id_k), intent(out) :: fid
      logical, intent(out) :: ok;  fid = 0;  ok = .false.
   end subroutine
   subroutine h5_close_file(fid);  integer(h5id_k), intent(in) :: fid;  end subroutine
   subroutine h5_group(loc, name, gid, ok)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      integer(h5id_k), intent(out) :: gid;  logical, intent(out) :: ok
      gid = 0;  ok = .false.
   end subroutine
   subroutine h5_group_open(loc, name, gid, ok)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      integer(h5id_k), intent(out) :: gid;  logical, intent(out) :: ok
      gid = 0;  ok = .false.
   end subroutine
   subroutine h5_group_close(gid);  integer(h5id_k), intent(in) :: gid;  end subroutine
   subroutine h5_unlink(loc, name)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
   end subroutine
   logical function h5_has(loc, name)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      h5_has = .false.
   end function
   subroutine h5_write_1d(loc, name, a, units, long_name)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), intent(in) :: a(:)
      character(len=*), intent(in), optional :: units, long_name
   end subroutine
   subroutine h5_write_2d(loc, name, a, units, long_name)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), intent(in) :: a(:,:)
      character(len=*), intent(in), optional :: units, long_name
   end subroutine
   subroutine h5_write_3d(loc, name, a, units, long_name)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), intent(in) :: a(:,:,:)
      character(len=*), intent(in), optional :: units, long_name
   end subroutine
   subroutine h5_write_6d(loc, name, a, units, long_name)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), intent(in) :: a(:,:,:,:,:,:)
      character(len=*), intent(in), optional :: units, long_name
   end subroutine
   subroutine h5_dims(loc, name, d, ok)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      integer(int64), allocatable, intent(out) :: d(:)
      logical, intent(out) :: ok
      ok = .false.
   end subroutine
   subroutine h5_read_1d(loc, name, a, ok, i0)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), allocatable, intent(out) :: a(:);  logical, intent(out) :: ok
      integer, intent(in), optional :: i0
      ok = .false.
   end subroutine
   subroutine h5_read_2d(loc, name, a, ok, i0)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), allocatable, intent(out) :: a(:,:);  logical, intent(out) :: ok
      integer, intent(in), optional :: i0
      ok = .false.
   end subroutine
   subroutine h5_read_3d(loc, name, a, ok, i0)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), allocatable, intent(out) :: a(:,:,:);  logical, intent(out) :: ok
      integer, intent(in), optional :: i0
      ok = .false.
   end subroutine
   subroutine h5_read_6d(loc, name, a, ok)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), allocatable, intent(out) :: a(:,:,:,:,:,:);  logical, intent(out) :: ok
      ok = .false.
   end subroutine
   subroutine h5_put_attr_d(loc, name, v)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), intent(in) :: v
   end subroutine
   subroutine h5_put_attr_i(loc, name, v)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      integer, intent(in) :: v
   end subroutine
   subroutine h5_put_attr_s(loc, name, v)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name, v
   end subroutine
   subroutine h5_get_attr_d(loc, name, v, ok)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      real(real64), intent(out) :: v;  logical, intent(out) :: ok
      v = 0.0_real64;  ok = .false.
   end subroutine
   subroutine h5_get_attr_i(loc, name, v, ok)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      integer, intent(out) :: v;  logical, intent(out) :: ok
      v = 0;  ok = .false.
   end subroutine
   subroutine h5_get_attr_s(loc, name, v, ok)
      integer(h5id_k), intent(in) :: loc;  character(len=*), intent(in) :: name
      character(len=*), intent(out) :: v;  logical, intent(out) :: ok
      v = '';  ok = .false.
   end subroutine
#endif

end module sedust_h5
