module sedust_product_mod
   ! Reads the SEDust optics products out of data/<model>/sedust_<model>.h5 -- the
   ! wavelength axis, the (lambda, a_eff) cross-section tables, and the
   ! size-integrated extinction curve.  sedust_h5.f90 knows HDF5; this knows
   ! where SEDust puts things inside one, and is the only place that does on
   ! the reading side (sed/calc_qtable.f90 and sed/calc_kext.f90 are the two on
   ! the writing side).
   !
   ! WHAT include_euv MEANS.  A file holds ONE wavelength axis, the widest the
   ! model has, and /grid records i_lyman, where that axis crosses the Lyman
   ! limit (lyman_index below defines exactly where).  Every routine here takes
   ! include_euv:
   !
   !   .true.   the whole axis, for a host that transports ionizing radiation
   !   .false.  lambda(i_lyman:) and the same slice of every wavelength-indexed
   !            array -- the non-ionizing product, which used to be a second
   !            file that had to be kept in step with the first
   !
   ! The slice is taken with an HDF5 hyperslab, so the narrow read moves only
   ! the rows it returns.
   !
   ! NOTHING HERE STOPS THE RUN.  Every entry point reports failure through ok,
   ! exactly as the text readers beside it do (kext_table.f90, q_component.f90),
   ! so a caller can fall back to the text products and say which it used.  A
   ! build without HDF5 gets ok = .false. from the layer below and needs no
   ! #ifdef of its own.
   use, intrinsic :: iso_fortran_env, only: real64
   use constants, only: wp
   use sedust_h5
   implicit none
   private

   ! The Lyman limit, 13.6 eV: where an interstellar radiation field stops, and
   ! the short-wavelength end of the non-ionizing view of every product.
   real(wp), parameter :: LAM_LYMAN_UM = 0.0912_wp
   ! A node this close to the limit in relative terms IS the limit.  The
   ! astrodust and DL07 axes carry 0.0912 um bit-exactly, so this changes
   ! nothing for them today; it keeps the rule from turning on the last bit of
   ! a product someone rewrites.
   real(wp), parameter :: LYMAN_TOL = 1.0e-9_wp

   public :: sedust_dir, sedust_h5_file
   public :: LAM_LYMAN_UM, lyman_index
   public :: read_sedust_grid
   public :: read_sedust_qtable
   public :: read_sedust_kext

contains

   pure integer function lyman_index(lam) result(i)
      ! Where an ascending wavelength axis crosses the Lyman limit: the LAST
      ! node at or below 0.0912 um, not the first at or above it.
      !
      ! The non-ionizing view of a product is lam(i:), and a host interpolates
      ! that view onto its own grid, so the view has to COVER the band the host
      ! transports.  A host whose transport floor IS the Lyman limit is then
      ! inside the table by construction, for every model alike.
      !
      ! Cutting at the first node at or above the limit instead leaves the
      ! limit OUTSIDE the view whenever the grid has no node exactly there.
      ! The astrodust and DL07 axes do have one -- both rules return the same
      ! index for them -- but the ZDA grid does not: its first node above the
      ! limit sits 1.2e-5 of itself inside it, just past dust_extinction's
      ! EDGE_TOL, so the same call that was served for two models was refused
      ! for the third.  One node of the model's own table is not extrapolation;
      ! widening the allowance instead would have weakened the guard for every
      ! model at every wavelength.
      real(wp), intent(in) :: lam(:)
      integer :: k
      i = 1
      do k = 1, size(lam)
         if (lam(k) > LAM_LYMAN_UM*(1.0_wp + LYMAN_TOL)) exit
         i = k
      end do
   end function lyman_index


   function sedust_dir(data_dir, model) result(d)
      ! Everything a model owns lives in ONE directory, data_dir/<model>/:
      ! its HDF5 product, its cross-section tables, its extinction curves and,
      ! where the model IS a set of files (ZDA), its definition.  What is shared
      ! between models -- the dielectric functions, the published reference
      ! tables -- stays beside those directories rather than inside one of them,
      ! because a D03 astrosilicate index is not DL07's any more than it is
      ! Zubko's.  A host then ships data_dir/<model>/ for what it runs.
      character(len=*), intent(in) :: data_dir, model
      character(len=512) :: d
      integer :: k
      k = len_trim(data_dir)
      if (k == 0) then
         d = trim(model)//'/'
      else if (data_dir(k:k) == '/') then
         d = trim(data_dir)//trim(model)//'/'
      else
         d = trim(data_dir)//'/'//trim(model)//'/'
      end if
   end function sedust_dir


   function sedust_h5_file(data_dir, model) result(path)
      ! The model's product, inside its own directory.
      character(len=*), intent(in) :: data_dir, model
      character(len=512) :: path
      path = trim(sedust_dir(data_dir, model))//'sedust_'//trim(model)//'.h5'
   end function sedust_h5_file


   subroutine read_sedust_grid(path, include_euv, lam, i_lyman, ok)
      ! The model's wavelength axis, cut as include_euv asks.  i_lyman comes
      ! back as the index it has in the FILE, which is what the qtable and kext
      ! readers below must be given to slice identically; on the returned
      ! narrow axis the corresponding index is 1.
      character(len=*),      intent(in)  :: path
      logical,               intent(in)  :: include_euv
      real(wp), allocatable, intent(out) :: lam(:)
      integer,               intent(out) :: i_lyman
      logical,               intent(out) :: ok

      integer(h5id_k) :: fid, gid
      logical :: got

      ok = .false.;  i_lyman = 1
      call h5_begin(got);  if (.not. got) return
      call h5_open_file(path, fid, got)
      if (.not. got) then;  call h5_end();  return;  end if
      call h5_group_open(fid, 'grid', gid, got)
      if (.not. got) then
         call h5_close_file(fid);  call h5_end();  return
      end if
      call h5_get_attr_i(gid, 'i_lyman', i_lyman, got)
      if (.not. got) i_lyman = 1
      call h5_read_1d(gid, 'lambda', lam, ok, i0 = lyman_offset(include_euv, i_lyman))
      call h5_group_close(gid)
      call h5_close_file(fid)
      call h5_end()
      if (.not. ok) then
         if (allocated(lam)) deallocate(lam)
         i_lyman = 1
      end if
   end subroutine read_sedust_grid


   subroutine read_sedust_qtable(path, comp, include_euv, aeff, Qext, Qabs, Qsca, &
                                 gpar, rho, ok, flag)
      ! One grain population's cross sections, (n_lam, n_a) on the cut axis.
      !
      ! A population that only absorbs -- the PAHs, whose cross sections are a
      ! prescription and not a Mie solution -- carries Q_abs alone in the file.
      ! Qext, Qsca and gpar then come back as ZERO rather than absent, which is
      ! what the size integral must use for it, and the group's `scatters`
      ! attribute records that this is meant rather than missing.
      character(len=*),      intent(in)  :: path, comp
      logical,               intent(in)  :: include_euv
      real(wp), allocatable, intent(out) :: aeff(:)
      real(wp), allocatable, intent(out) :: Qext(:,:), Qabs(:,:), Qsca(:,:), gpar(:,:)
      ! Solid density of the material [g/cm^3]; 0 when the group states none.
      real(wp),              intent(out) :: rho
      logical,               intent(out) :: ok
      ! The astrodust table's regime flag: 0 T-matrix, 10 Rayleigh dipole,
      ! 20 geometric optics.  Left unallocated for a group that has none.
      integer, allocatable, optional, intent(out) :: flag(:,:)

      integer(h5id_k) :: fid, gid
      real(wp), allocatable :: rflag(:,:)
      logical :: got
      integer :: i_lyman, i0

      ok = .false.;  rho = 0.0_wp
      call h5_begin(got);  if (.not. got) return
      call h5_open_file(path, fid, got)
      if (.not. got) then;  call h5_end();  return;  end if

      call read_lyman_index(fid, i_lyman)
      i0 = lyman_offset(include_euv, i_lyman)

      call h5_group_open(fid, 'qtable/'//trim(comp), gid, got)
      if (.not. got) then
         call h5_close_file(fid);  call h5_end();  return
      end if

      call h5_read_1d(gid, 'a_eff', aeff, ok)
      if (.not. ok) then
         call h5_group_close(gid);  call h5_close_file(fid);  call h5_end();  return
      end if
      call h5_read_2d(gid, 'Q_abs', Qabs, ok, i0 = i0)
      if (.not. ok) then
         deallocate(aeff)
         call h5_group_close(gid);  call h5_close_file(fid);  call h5_end();  return
      end if

      if (h5_has(gid, 'Q_sca')) then
         call h5_read_2d(gid, 'Q_ext', Qext, got, i0 = i0);  ok = ok .and. got
         call h5_read_2d(gid, 'Q_sca', Qsca, got, i0 = i0);  ok = ok .and. got
         call h5_read_2d(gid, 'g',     gpar, got, i0 = i0);  ok = ok .and. got
      else
         allocate(Qext(size(Qabs,1), size(Qabs,2)))
         allocate(Qsca(size(Qabs,1), size(Qabs,2)))
         allocate(gpar(size(Qabs,1), size(Qabs,2)))
         Qext = Qabs;  Qsca = 0.0_wp;  gpar = 0.0_wp
      end if

      if (present(flag)) then
         if (h5_has(gid, 'regime')) then
            call h5_read_2d(gid, 'regime', rflag, got, i0 = i0)
            if (got) then
               allocate(flag(size(rflag,1), size(rflag,2)))
               flag = nint(rflag)
            end if
            if (allocated(rflag)) deallocate(rflag)
         end if
      end if

      call h5_get_attr_d(gid, 'rho_bulk_g_cm3', rho, got)
      if (.not. got) rho = 0.0_wp

      call h5_group_close(gid)
      call h5_close_file(fid)
      call h5_end()
      if (.not. ok) call drop_qtable()

   contains

      subroutine drop_qtable()
         if (allocated(aeff)) deallocate(aeff)
         if (allocated(Qabs)) deallocate(Qabs)
         if (allocated(Qext)) deallocate(Qext)
         if (allocated(Qsca)) deallocate(Qsca)
         if (allocated(gpar)) deallocate(gpar)
         if (present(flag)) then
            if (allocated(flag)) deallocate(flag)
         end if
         rho = 0.0_wp
      end subroutine drop_qtable

   end subroutine read_sedust_qtable


   subroutine read_sedust_kext(path, include_euv, n, lam, albedo, gbar, &
                               Cext, Cabs, Csca, ok, Mdust_H, group)
      ! The size-integrated extinction curve, in the argument order
      ! load_kext_table uses, so that a caller can take either source without
      ! reshaping anything.  albedo is the file's own C_sca/C_ext.
      !
      ! group names the curve inside the file, '/kext' by default.  A model
      ! that stores more than one set of optics stores the matching curve
      ! beside each -- zubko has /kext for the distributed tables and
      ! /kext_mie_d03 for the recomputation -- so that the curve a host is
      ! served is the size integral of the very optics its model was built on.
      character(len=*),      intent(in)  :: path
      logical,               intent(in)  :: include_euv
      integer,               intent(out) :: n
      real(wp), allocatable, intent(out) :: lam(:), albedo(:), gbar(:)
      real(wp), allocatable, intent(out) :: Cext(:), Cabs(:), Csca(:)
      logical,               intent(out) :: ok
      ! Dust mass per H [g/H] the K_abs column is normalized by; 0 when the
      ! group states none.
      real(wp), optional,    intent(out) :: Mdust_H
      character(len=*), optional, intent(in) :: group

      integer(h5id_k) :: fid, gid
      logical :: got
      integer :: i_lyman, i0
      character(len=64) :: grp

      grp = 'kext';  if (present(group)) grp = group
      n = 0;  ok = .false.
      if (present(Mdust_H)) Mdust_H = 0.0_wp
      call h5_begin(got);  if (.not. got) return
      call h5_open_file(path, fid, got)
      if (.not. got) then;  call h5_end();  return;  end if

      call read_lyman_index(fid, i_lyman)
      i0 = lyman_offset(include_euv, i_lyman)

      call h5_group_open(fid, trim(grp), gid, got)
      if (.not. got) then
         ! The file exists but carries no curve: calc_qtable.x replaces the file
         ! and calc_kext.x has not run since.  Absence must be tellable from a
         ! computed zero, so this is a failure, not an empty result.
         call h5_close_file(fid);  call h5_end();  return
      end if

      call h5_read_1d(fid, 'grid/lambda', lam, ok, i0 = i0)
      if (ok) call h5_read_1d(gid, 'albedo', albedo, ok, i0 = i0)
      if (ok) call h5_read_1d(gid, 'g',      gbar,   ok, i0 = i0)
      if (ok) call h5_read_1d(gid, 'C_ext',  Cext,   ok, i0 = i0)
      if (ok) call h5_read_1d(gid, 'C_abs',  Cabs,   ok, i0 = i0)
      if (ok) call h5_read_1d(gid, 'C_sca',  Csca,   ok, i0 = i0)
      if (ok .and. present(Mdust_H)) then
         call h5_get_attr_d(gid, 'M_dust_per_H', Mdust_H, got)
         if (.not. got) Mdust_H = 0.0_wp
      end if

      call h5_group_close(gid)
      call h5_close_file(fid)
      call h5_end()

      if (ok) then
         n = size(lam)
         ! The same two conditions load_kext_table refuses on, so that neither
         ! source can hand the interpolation a grid it cannot bracket in.
         if (n < 2) ok = .false.
         if (ok) ok = strictly_ascending(lam)
      end if
      if (.not. ok) then
         if (allocated(lam))    deallocate(lam)
         if (allocated(albedo)) deallocate(albedo)
         if (allocated(gbar))   deallocate(gbar)
         if (allocated(Cext))   deallocate(Cext)
         if (allocated(Cabs))   deallocate(Cabs)
         if (allocated(Csca))   deallocate(Csca)
         n = 0
      end if
   end subroutine read_sedust_kext


   ! ===================================================================
   subroutine read_lyman_index(fid, i_lyman)
      ! /grid's i_lyman, defaulting to 1 -- the whole axis -- for a file that
      ! states none, which is the reading under which a narrow request returns
      ! everything rather than silently dropping the short end.
      integer(h5id_k), intent(in)  :: fid
      integer,         intent(out) :: i_lyman
      integer(h5id_k) :: gid
      logical :: got
      i_lyman = 1
      call h5_group_open(fid, 'grid', gid, got);  if (.not. got) return
      call h5_get_attr_i(gid, 'i_lyman', i_lyman, got)
      if (.not. got) i_lyman = 1
      call h5_group_close(gid)
   end subroutine read_lyman_index


   pure integer function lyman_offset(include_euv, i_lyman) result(i0)
      logical, intent(in) :: include_euv
      integer, intent(in) :: i_lyman
      i0 = 1
      if (.not. include_euv) i0 = max(1, i_lyman)
   end function lyman_offset


   pure logical function strictly_ascending(x)
      real(wp), intent(in) :: x(:)
      integer :: i
      strictly_ascending = .true.
      do i = 2, size(x)
         if (x(i) <= x(i-1)) then
            strictly_ascending = .false.;  return
         end if
      end do
   end function strictly_ascending

end module sedust_product_mod
