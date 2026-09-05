program calc_polarized_optics
   ! Writes the /polarized groups of data/astrodust/sedust_astrodust.h5: the
   ! orientation-resolved cross sections of the DH21 astrodust spheroid, and
   ! the random-orientation and aligned scattering matrices.
   !
   !   ./calc_polarized_optics.x [qjori | scatmat | all]        (default: all)
   !
   ! These are the products of the polarized branch, and they belong in the
   ! astrodust model's own file beside its scalar optics: one file per model is
   ! what a host ships.  This program APPENDS, exactly as calc_kext.x does --
   ! calc_qtable.x replaces the file, so the running order is
   !
   !   ./calc_qtable.x astrodust          (lays down /grid and /qtable)
   !   ./calc_kext.x astrodust euv        (adds /kext)
   !   ./calc_polarized_optics.x          (adds /polarized)
   !
   ! WHY /polarized CARRIES ITS OWN WAVELENGTH AXIS.  The scalar products run
   ! into the ionizing band; the orientation-resolved table does not.  It stops
   ! at the Lyman limit unless the EUV companion table has been computed, and
   ! below its first node build_Cpol leaves the dichroic extinction at exactly
   ! zero -- an omission, not physics, since a b/a = 1.4 spheroid has a nonzero
   ! dichroic extinction there.  So /polarized is filed against
   ! /polarized/lambda rather than /grid/lambda, and records covers_euv and
   ! pol_valid_from, so that a reader can tell a band that was not computed
   ! from one that came out zero.  It is NOT sliced by /grid's i_lyman.
   use, intrinsic :: iso_fortran_env, only: real64
   use constants,          only: wp
   use sedust_h5
   use q_table_jori_mod,   only: read_jori_stream, A_ALIGN, ALPHA_ALIGN, FMAX_ALIGN
   use scatmat_aligned_mod, only: load_scatmat_aligned, free_scatmat_aligned, &
                                 scm_loaded, scm_nband, scm_nti, scm_nts, scm_nphi, &
                                 scm_lambda, scm_theta_i, scm_theta_s, &
                                 scm_phi, scm_theta_ran, &
                                 scm_cext_al, scm_cpol_al, scm_cbir_al, scm_csca_al, &
                                 scm_csca_pol_al, scm_cext_tot, scm_csca_tot, &
                                 scm_cext_ref, scm_csca_ref, scm_F_tot, scm_F_ref, &
                                 scm_Z, scm_profile_name, scm_fmax, scm_a_align, &
                                 scm_alpha
   use sed_run_options, only: run_options_t, declare_run_options, &
                              read_run_subject, read_run_option
   implicit none

   character(len=*), parameter :: H5FILE  = '../data/astrodust/sedust_astrodust.h5'
   character(len=*), parameter :: F_QJ    = '../data/astrodust/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: F_WAVE  = '../data/astrodust/DH21_wave'
   character(len=*), parameter :: F_AEFF  = '../data/astrodust/DH21_aeff'
   character(len=*), parameter :: F_SCAT  = &
        '../data/astrodust/scatmat_aligned_astrodust_P0.20_Fe0.00_1.400.dat.gz'
   ! The DH21 table's own grid lengths, as q_table_jori_mod states them.
   integer, parameter :: NA_J = 169, NW_J = 1129
   real(wp), parameter :: LAM_LYMAN = 0.0912_wp

   character(len=32) :: which
   integer :: narg, iarg
   logical :: taken
   type(run_options_t) :: o
   character(len=64)   :: arg

   ! Which polarized product to write is the only axis here.  Both products are
   ! read from the DH21 tables on those tables' own grids, so there is no
   ! wavelength-grid axis, no field and no solver; a word naming one of those
   ! is refused by name.
   call declare_run_options(o, program='calc_polarized_optics', &
        subjects=[character(len=16):: 'qjori', 'scatmat', 'all'])
   narg = command_argument_count()
   which = 'all'
   if (narg >= 1) then
      call get_command_argument(1, which)
      call read_run_subject(o, trim(which), taken)
      if (.not. taken) then
         write(*,'(a,a)') ' calc_polarized_optics: unknown product ', trim(which)
         write(*,'(a)')   ' usage: ./calc_polarized_optics.x [qjori | scatmat | all]'
         stop 1
      end if
   end if
   do iarg = 2, narg
      call get_command_argument(iarg, arg)
      call read_run_option(trim(arg), o, taken)
      if (.not. taken) then
         write(*,'(a,a)') ' calc_polarized_optics: unknown argument ', trim(arg)
         write(*,'(a)')   ' usage: ./calc_polarized_optics.x [qjori | scatmat | all]'
         stop 1
      end if
   end do

   if (.not. sedust_has_hdf5) then
      write(*,'(a)') ' calc_polarized_optics: built without HDF5;'// &
                     ' there is nothing to write'
      stop 1
   end if

   select case (trim(which))
   case ('qjori');    call write_qjori()
   case ('scatmat');  call write_scatmat()
   case ('all');      call write_qjori();  call write_scatmat()
   end select

contains

   ! ===================================================================
   subroutine write_qjori()
      ! The orientation-resolved table, all three orientations kept:
      !   jori = 1  k || a          (a = the spheroid symmetry axis)
      !   jori = 2  k perp a, E || a
      !   jori = 3  k perp a, E perp a
      ! The reduced combinations the solver uses -- the dichroic
      ! 0.5*(Q(3) - Q(2)) and the random-orientation average -- are NOT stored:
      ! they are one subtraction away from these, and storing a derived
      ! quantity beside the thing it is derived from is how the two drift
      ! apart.
      real(wp), allocatable :: lam(:), aeff(:)
      real(wp), allocatable :: qe(:,:,:), qa(:,:,:), qs(:,:,:), qr(:,:,:)
      logical :: bir, ok
      character(len=256) :: msg
      integer(h5id_k) :: fid, gid

      call read_jori_stream(F_QJ, F_WAVE, F_AEFF, NA_J, NW_J, &
                            lam, aeff, qe, qa, qs, qr, bir, ok, msg)
      if (.not. ok) then
         write(*,'(a,a)') ' calc_polarized_optics: ', trim(msg)
         stop 1
      end if
      write(*,'(a,i0,a,i0,a,l1)') ' qjori: ', size(lam), ' wavelengths x ', &
         size(aeff), ' radii x 3 orientations, Q_re present: ', bir

      call open_polarized(fid, ok, lam)
      if (.not. ok) stop 1
      if (h5_has(fid, 'polarized/qjori')) call h5_unlink(fid, 'polarized/qjori')
      call h5_group(fid, 'polarized/qjori', gid, ok)
      if (.not. ok) then
         write(*,'(a)') ' calc_polarized_optics: cannot create /polarized/qjori'
         call h5_close_file(fid);  call h5_end();  stop 1
      end if

      call h5_write_3d_jori(gid, 'Q_ext', qe, single=.true.)
      call h5_write_3d_jori(gid, 'Q_abs', qa, single=.true.)
      call h5_write_3d_jori(gid, 'Q_sca', qs, single=.true.)
      if (bir) call h5_write_3d_jori(gid, 'Q_re', qr, single=.true.)
      call h5_write_1d(gid, 'a_eff', aeff, units='um', &
                       long_name='grain effective radius')
      call h5_put_attr_s(gid, 'jori_convention', &
           '1 = k || a; 2 = k perp a, E || a; 3 = k perp a, E perp a '// &
           '(a is the spheroid symmetry axis)')
      call h5_put_attr_s(gid, 'method', &
           'fixed-orientation T-matrix on the DH21 astrodust oblate spheroid, b/a = 1.400')
      call h5_put_attr_s(gid, 'source', F_QJ)
      ! The alignment profile the SOLVER applies to these cross sections.  It
      ! is not the one the aligned scattering matrix was integrated under --
      ! that one is recorded on its own group -- and the two can differ.
      call h5_put_attr_d(gid, 'falign_f_max',   FMAX_ALIGN)
      call h5_put_attr_d(gid, 'falign_a_align', A_ALIGN)
      call h5_put_attr_d(gid, 'falign_alpha',   ALPHA_ALIGN)
      call h5_group_close(gid)
      call h5_close_file(fid)
      call h5_end()
      write(*,'(a,a)') ' wrote /polarized/qjori into ', H5FILE
      deallocate(lam, aeff, qe, qa, qs)
      if (allocated(qr)) deallocate(qr)
   end subroutine write_qjori


   ! ===================================================================
   subroutine write_scatmat()
      ! The scattering matrices, in the two forms the table carries: the
      ! random-orientation F on its own scattering-angle grid, and the aligned
      ! phase matrix Z on the (theta_i, theta_s, phi) grid the polarized
      ! transfer needs.  They are separate groups because they are separate
      ! products -- F is what an unaligned population scatters, Z is what an
      ! aligned one does at the reference alignment eta = 1 -- and a host uses
      ! one, the other, or a blend of the two.
      integer(h5id_k) :: fid, gid
      logical :: ok
      integer :: stat

      call load_scatmat_aligned(F_SCAT, stat)
      if (stat /= 0 .or. .not. scm_loaded) then
         write(*,'(a,a,a,i0)') ' calc_polarized_optics: cannot read ', F_SCAT, ', status = ', stat
         stop 1
      end if
      write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' scatmat: ', scm_nband, ' bands, Z is (', &
         scm_nti, ', ', scm_nts, ', ', scm_nphi, ', 4, 4, nband)'

      call open_polarized(fid, ok)
      if (.not. ok) stop 1

      ! ---- random orientation ------------------------------------------
      if (h5_has(fid, 'polarized/scatmat_random')) &
         call h5_unlink(fid, 'polarized/scatmat_random')
      if (h5_has(fid, 'polarized/scatmat_aligned')) &
         call h5_unlink(fid, 'polarized/scatmat_aligned')
      call h5_group(fid, 'polarized/scatmat_random', gid, ok)
      if (.not. ok) then;  call h5_close_file(fid);  call h5_end();  stop 1;  end if
      call h5_write_1d(gid, 'band_lambda', scm_lambda, units='um', &
                       long_name='band centre wavelength')
      call h5_write_1d(gid, 'theta', scm_theta_ran, units='deg', &
                       long_name='scattering angle of the F matrix')
      call h5_write_3d(gid, 'F_tot', scm_F_tot, units='1', &
           long_name='random-orientation scattering matrix, total population '// &
                     '(6 elements: 11, 22, 33, 44, 12, 34)', single=.true.)
      call h5_write_3d(gid, 'F_ref', scm_F_ref, units='1', &
           long_name='the same for the reference (aligned-grain) population alone', &
           single=.true.)
      call h5_write_1d(gid, 'C_sca_tot', scm_csca_tot, units='um^2/H', &
                       long_name='scattering cross section per H, total population', &
                       single=.true.)
      call h5_write_1d(gid, 'C_sca_ref', scm_csca_ref, units='um^2/H', &
                       long_name='the same for the reference population', single=.true.)
      call h5_write_1d(gid, 'C_ext_tot', scm_cext_tot, units='um^2/H', &
                       long_name='extinction cross section per H, total population', &
                       single=.true.)
      call h5_write_1d(gid, 'C_ext_ref', scm_cext_ref, units='um^2/H', &
                       long_name='the same for the reference population', single=.true.)
      call h5_put_attr_s(gid, 'F_normalization', &
           'alpha1-normalized: (1/2) INT F11 dcos(theta) = 1, i.e. INT F11 dOmega = 4 pi. '// &
           'The absolute differential matrix [um^2 sr^-1 per H] is C_sca * F / (4 pi).')
      call h5_put_attr_s(gid, 'source', F_SCAT)
      call h5_group_close(gid)

      ! ---- aligned ------------------------------------------------------
      call h5_group(fid, 'polarized/scatmat_aligned', gid, ok)
      if (.not. ok) then;  call h5_close_file(fid);  call h5_end();  stop 1;  end if
      call h5_write_1d(gid, 'band_lambda', scm_lambda, units='um', &
                       long_name='band centre wavelength')
      call h5_write_1d(gid, 'theta_i', scm_theta_i, units='deg', &
                       long_name='angle between the incident direction and the alignment axis')
      call h5_write_1d(gid, 'theta_s', scm_theta_s, units='deg', &
                       long_name='scattering polar angle')
      call h5_write_1d(gid, 'phi', scm_phi, units='deg', &
                       long_name='scattering azimuth')
      call h5_write_6d(gid, 'Z', scm_Z, units='um^2 sr^-1 per H', &
           long_name='aligned phase matrix at the reference alignment eta = 1', &
           single=.true.)
      call h5_write_2d(gid, 'C_ext_al',     scm_cext_al,     units='um^2/H', &
                       long_name='extinction cross section per H of the aligned population', &
                       single=.true.)
      call h5_write_2d(gid, 'C_pol_al',     scm_cpol_al,     units='um^2/H', &
                       long_name='dichroic (polarized) extinction cross section per H', &
                       single=.true.)
      call h5_write_2d(gid, 'C_bir_al',     scm_cbir_al,     units='um^2/H', &
                       long_name='birefringent (circular) extinction cross section per H', &
                       single=.true.)
      call h5_write_2d(gid, 'C_sca_al',     scm_csca_al,     units='um^2/H', &
                       long_name='scattering cross section per H, by grid closure over Z', &
                       single=.true.)
      call h5_write_2d(gid, 'C_sca_pol_al', scm_csca_pol_al, units='um^2/H', &
                       long_name='polarized part of the aligned scattering cross section', &
                       single=.true.)
      ! The alignment profile this table was INTEGRATED under.  A model whose
      ! own profile differs is using a stale table, which is the check
      ! alignment_matches_scatmat makes at run time and which these attributes
      ! let a reader make at open time.
      call h5_put_attr_s(gid, 'profile_name', trim(scm_profile_name))
      call h5_put_attr_d(gid, 'f_max',   scm_fmax)
      call h5_put_attr_d(gid, 'a_align', scm_a_align)
      call h5_put_attr_d(gid, 'alpha',   scm_alpha)
      call h5_put_attr_s(gid, 'eta_meaning', &
           'Z and the aligned cross sections are at eta = 1; a host scales them '// &
           'by its own alignment strength.')
      call h5_put_attr_s(gid, 'source', F_SCAT)
      call h5_group_close(gid)

      call h5_close_file(fid)
      call h5_end()
      write(*,'(a,a)') ' wrote /polarized/scatmat_{random,aligned} into ', H5FILE
      call free_scatmat_aligned()
   end subroutine write_scatmat


   ! ===================================================================
   subroutine open_polarized(fid, ok, lam)
      ! Open the model's file for appending and make sure /polarized exists.
      !
      ! The axis is laid down by the FIRST writer that has one, which is the
      ! orientation-resolved table: /polarized/lambda is its wavelength grid.
      ! The scattering matrices are banded and carry their own band_lambda on
      ! each group, so they pass no lam and never define the axis -- running
      ! `calc_polarized_optics.x scatmat` on a fresh file must not leave /polarized
      ! claiming a five-point wavelength axis.
      integer(h5id_k),    intent(out) :: fid
      logical,            intent(out) :: ok
      real(real64), optional, intent(in) :: lam(:)
      integer(h5id_k) :: pid

      call h5_begin(ok);  if (.not. ok) return
      call h5_open_rw(H5FILE, fid, ok)
      if (.not. ok) then
         write(*,'(a,a)') ' calc_polarized_optics: cannot open ', H5FILE
         write(*,'(a)')   '   write it first with  ./calc_qtable.x astrodust'
         call h5_end();  return
      end if

      if (.not. h5_has(fid, 'polarized')) then
         call h5_group(fid, 'polarized', pid, ok)
         if (.not. ok) then
            call h5_close_file(fid);  call h5_end();  return
         end if
         if (present(lam)) then
            call h5_write_1d(pid, 'lambda', lam, units='um', &
                 long_name='wavelength axis of the orientation-resolved product')
            ! covers_euv = 0 says the ionizing band is ABSENT from this product,
            ! not that it is zero there.  pol_valid_from is where the table
            ! actually starts.
            call h5_put_attr_i(pid, 'covers_euv', merge(1, 0, lam(1) < LAM_LYMAN))
            call h5_put_attr_d(pid, 'pol_valid_from', lam(1))
            call h5_put_attr_s(pid, 'covers_euv_meaning', &
                 '0 = the ionizing band is absent from this product; a dichroic '// &
                 'extinction of zero below pol_valid_from is an omission, not physics')
         end if
         call h5_put_attr_s(pid, 'generator', 'SEDust sed/calc_polarized_optics.x')
         call h5_group_close(pid)
      end if
   end subroutine open_polarized


   subroutine h5_write_3d_jori(gid, name, q, single)
      ! (n_lam, n_a, 3) in Fortran, which h5py reads as (3, n_a, n_lam) --
      ! the shape section 2 of the migration note states.
      integer(h5id_k),  intent(in) :: gid
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: q(:,:,:)
      logical,          intent(in), optional :: single
      call h5_write_3d(gid, name, q, units='1', &
           long_name='efficiency factor, one block per orientation (see jori_convention)', &
           single=single)
   end subroutine h5_write_3d_jori


end program calc_polarized_optics
