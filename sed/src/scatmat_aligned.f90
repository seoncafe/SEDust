module scatmat_aligned_mod
   ! Size-integrated fixed-orientation (aligned) scattering optics of the
   ! DH21 astrodust spheroid, read from the ASCII product written by
   ! tmatrix/driver/run_scatmat_aligned.x and served to a polarized
   ! radiative-transfer host along a photon path.
   !
   ! The interface follows the plan (docs/tmatrix_aligned_scattering_plan.md,
   ! Section 4): MoCafe supplies the alignment assumptions and the local field
   ! geometry, SEDust returns the matrices, MoCafe carries out the polarized
   ! transfer (Peest et al. 2017/2023 formalism). This module owns only the
   ! material quantities in matrix form.
   !
   ! TWO STRICTLY SEPARATED LAYERS
   !
   !   INITIALIZATION (serial, called once before the photon loop):
   !     load_scatmat_aligned -- parse the table into module storage, form the
   !     cos(theta_s) abscissae, record the alignment profile. NOT thread-safe.
   !
   !   PATH QUERIES (pure reads, no allocation, no I/O, no state mutation):
   !     scatmat_band, extinction_matrix_aligned, mueller_matrix_aligned,
   !     mueller_matrix_random, scattering_cross_sections. Each writes only its
   !     own output arguments and reads the loaded arrays, so it is SAFE to call
   !     concurrently from OpenMP photon threads. The per-cell dependence enters
   !     only through the scalars the host already holds: theta_i (from one dot
   !     product k-hat . B-hat) and the alignment scale eta.
   !
   ! The loaded grids and arrays are exposed read-only (public, protected) so a
   ! host can take them once at startup and inline its own lookups without
   ! calling back into the library during transport.
   !
   ! ETA CONTRACT (from the table header). For a cell alignment scale eta the
   ! aligned optics scale linearly: Z_al,cell = eta Z_al, K_al,cell = eta K_al,
   ! and the unaligned remainder scattering matrix is
   !   F_unal = Csca_tot F_tot - eta Csca_ref F_ref   (absolute units),
   ! the unaligned extinction adding the isotropic Cext_tot - Cext_ref
   ! (its Cpol = Cbir = 0). The linearity in f_align is exact.
   !
   ! STOKES BASIS. Mishchenko meridional (v,h) = (theta-hat, phi-hat) of each
   ! propagation direction in the grain frame (z = alignment axis), Q = Iv - Ih.
   !
   ! MEMORY. Dominated by scm_Z: nti*nts*nphi*16*nband doubles. For the
   ! production grid (nti=19, nts=181, nphi=37, nband=5) that is
   ! 19*181*37*16*5*8 bytes ~ 81 MB. The F and K arrays add < 0.1 MB.

   use, intrinsic :: iso_fortran_env, only: error_unit, int64
   use constants, only: wp, deg2rad
   implicit none
   private

   ! Initialization layer
   public :: load_scatmat_aligned, free_scatmat_aligned
   ! Alignment-consistency guard, called by dust_set_alignment*
   public :: alignment_matches_scatmat
   ! Query layer
   public :: scatmat_band, extinction_matrix_aligned, mueller_matrix_aligned, &
             mueller_matrix_random, scattering_cross_sections

   ! ---- read-only storage exposure (public, protected) ------------------
   logical,        protected, public :: scm_loaded = .false.
   integer,        protected, public :: scm_nband = 0
   ! The aligned phase matrix Z is on (theta_i, theta_s, phi); the random-
   ! orientation matrices F are on their own, finer, scattering-angle grid
   ! (Theta, 1-degree). scm_nts (Z theta_s) and scm_ntheta (F Theta) are NOT the
   ! same grid -- they coincide only when the table was generated with a Z
   ! theta_s step of 1 degree.
   integer,        protected, public :: scm_nti = 0, scm_nts = 0, scm_nphi = 0
   integer,        protected, public :: scm_ntheta = 0
   real(wp), allocatable, protected, public :: scm_lambda(:)      ! (nband) [um]
   real(wp), allocatable, protected, public :: scm_theta_i(:)     ! (nti)   [deg] Z incidence
   real(wp), allocatable, protected, public :: scm_theta_s(:)     ! (nts)   [deg] Z scattering
   real(wp), allocatable, protected, public :: scm_phi(:)         ! (nphi)  [deg] Z azimuth
   real(wp), allocatable, protected, public :: scm_theta_ran(:)   ! (ntheta)[deg] F scattering
   real(wp), allocatable, protected, public :: scm_cos_theta_s(:) ! (nts)   cos abscissae
   ! K block on the theta_i grid [um^2/H]
   real(wp), allocatable, protected, public :: scm_cext_al(:,:)   ! (nti, nband)
   real(wp), allocatable, protected, public :: scm_cpol_al(:,:)   ! (nti, nband)
   real(wp), allocatable, protected, public :: scm_cbir_al(:,:)   ! (nti, nband)
   real(wp), allocatable, protected, public :: scm_csca_al(:,:)   ! (nti, nband) grid closure
   ! Per-band scalars [um^2/H]
   real(wp), allocatable, protected, public :: scm_cext_tot(:), scm_csca_tot(:)  ! (nband)
   real(wp), allocatable, protected, public :: scm_cext_ref(:), scm_csca_ref(:)  ! (nband)
   ! Random-orientation matrices (alpha1-normalized: (1/2) INT F11 dcos = 1),
   ! six elements (11, 22, 33, 44, 12, 34), on the F Theta grid. Restore absolute
   ! um^2 sr^-1 per H by multiplying F_tot by scm_csca_tot and F_ref by scm_csca_ref.
   real(wp), allocatable, protected, public :: scm_F_tot(:,:,:)   ! (ntheta, 6, nband)
   real(wp), allocatable, protected, public :: scm_F_ref(:,:,:)   ! (ntheta, 6, nband)
   ! Aligned phase matrix, um^2 sr^-1 per H at the reference alignment (eta=1).
   real(wp), allocatable, protected, public :: scm_Z(:,:,:,:,:,:) ! (nti,nts,nphi,4,4,nband)
   ! Alignment profile the table was integrated under.
   character(len=64), protected, public :: scm_profile_name = ''
   real(wp),          protected, public :: scm_fmax = 0.0_wp
   real(wp),          protected, public :: scm_a_align = 0.0_wp   ! [um]
   real(wp),          protected, public :: scm_alpha = 0.0_wp
   ! Set .true. by alignment_matches_scatmat when the model's alignment profile
   ! departs from the one this table was integrated under.
   logical,        protected, public :: scm_profile_mismatch = .false.
   ! Bytes actually allocated by the load (measured, not estimated).
   integer(int64), protected, public :: scm_bytes = 0_int64

   ! Fractional tolerance for scatmat_band's "exact" verdict and for the
   ! alignment-profile comparison in alignment_matches_scatmat.
   real(wp), parameter :: BAND_TOL    = 1.0e-3_wp
   real(wp), parameter :: PROFILE_TOL = 1.0e-6_wp

contains

   ! ====================================================================
   ! Initialization layer
   ! ====================================================================

   subroutine load_scatmat_aligned(path, status)
      ! Parse the aligned-scattering ASCII table into module storage.
      !
      ! status (0 = success); errors follow the dust_lib convention -- reported
      ! through the argument, never a STOP inside the library:
      !   1  cannot open the file
      !   2  structure error (no bands, or grid sizes not self-consistent)
      !   3  read/parse error in a data block
      !   4  a grid axis differs between bands
      character(len=*), intent(in)  :: path
      integer,          intent(out) :: status

      integer :: u, ios
      logical :: gz, sub_ok
      character(len=512) :: read_path

      status = 0
      call free_scatmat_aligned()

      ! A compressed table is expanded to a scratch copy, read, then removed;
      ! a path that does not end in .gz is opened directly.
      gz = .false.
      if (len_trim(path) > 3) gz = (path(len_trim(path)-2:len_trim(path)) == '.gz')
      if (gz) then
         read_path = 'scatmat_aligned_scratch.dat'
         call gunzip_to(path, trim(read_path), sub_ok)
         if (.not. sub_ok) then
            call discard_scratch(.true., trim(read_path))
            status = 1;  return
         end if
      else
         read_path = path
      end if

      open(newunit=u, file=trim(read_path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call discard_scratch(gz, trim(read_path))
         status = 1;  return
      end if

      ! ---- pass 1: structure (band count and the three grid sizes) -------
      call scan_structure(u, status)
      if (status /= 0) then
         close(u);  call discard_scratch(gz, trim(read_path));  call free_scatmat_aligned()
         return
      end if

      call allocate_storage()

      ! ---- pass 2: fill ---------------------------------------------------
      rewind(u)
      call read_bands(u, status)
      close(u)
      call discard_scratch(gz, trim(read_path))
      if (status /= 0) then
         call free_scatmat_aligned()
         return
      end if

      ! Precompute the scattering-angle cosine abscissae so path queries do no
      ! trigonometry beyond the linear interpolation weights.
      scm_cos_theta_s = cos(scm_theta_s * deg2rad)

      scm_loaded = .true.
      scm_bytes  = storage_bytes()
   end subroutine load_scatmat_aligned


   subroutine scan_structure(u, status)
      ! One string-only scan to fix all the grid sizes:
      !   scm_nband  = "# lambda" band headers over the whole file
      !   scm_nti    = K data rows of the first band (theta_i grid)
      !   scm_ntheta = F data rows of the first band (F Theta grid)
      !   scm_nphi   = leading Z rows sharing the first (theta_i, theta_s)
      !   scm_nts    = (Z rows at the first theta_i) / scm_nphi
      ! The Z block is nested theta_i (outer) -> theta_s -> phi (inner), so the
      ! azimuth cycles first; that fixes nphi and then nts without assuming the
      ! Z theta_s grid equals the F Theta grid.
      integer, intent(in)  :: u
      integer, intent(out) :: status
      character(len=512) :: line, s, content
      integer :: ios, mode, kti, fts, zc, rows_ti, jos
      real(wp) :: c1, c2, ti0, ts0
      logical  :: first_z, counting_phi, counting_ti

      status = 0
      scm_nband = 0;  kti = 0;  fts = 0;  zc = 0;  mode = 0
      scm_nphi = 0;  rows_ti = 0
      ti0 = 0.0_wp;  ts0 = 0.0_wp
      first_z = .true.;  counting_phi = .true.;  counting_ti = .true.
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         s = adjustl(line)
         if (len_trim(s) == 0) cycle
         if (s(1:1) == '#') then
            mode = 0                       ! any comment ends a data block
            content = adjustl(s(2:))
            if (starts_with(content, 'lambda')) scm_nband = scm_nband + 1
            if (scm_nband == 1) then
               if (starts_with(content, 'K block')) mode = 1
               if (starts_with(content, 'F block')) mode = 2
               if (starts_with(content, 'Z block')) mode = 3
            end if
         else
            select case (mode)
            case (1);  kti = kti + 1
            case (2);  fts = fts + 1
            case (3)
               zc = zc + 1
               read(s, *, iostat=jos) c1, c2
               if (jos /= 0) then;  status = 3;  return;  end if
               if (first_z) then
                  ti0 = c1;  ts0 = c2;  first_z = .false.
               end if
               if (counting_phi) then
                  if (c1 == ti0 .and. c2 == ts0) then
                     scm_nphi = scm_nphi + 1
                  else
                     counting_phi = .false.
                  end if
               end if
               if (counting_ti) then
                  if (c1 == ti0) then
                     rows_ti = rows_ti + 1
                  else
                     counting_ti = .false.
                  end if
               end if
            end select
         end if
      end do

      scm_nti = kti;  scm_ntheta = fts
      if (scm_nband < 1 .or. scm_nti < 2 .or. scm_ntheta < 2 .or. scm_nphi < 2) then
         status = 2;  return
      end if
      if (mod(rows_ti, scm_nphi) /= 0) then
         status = 2;  return
      end if
      scm_nts = rows_ti / scm_nphi
      if (scm_nts < 2 .or. zc /= scm_nti*scm_nts*scm_nphi) status = 2
   end subroutine scan_structure


   subroutine read_bands(u, status)
      ! Second pass: parse the alignment profile, the per-band scalars, and the
      ! K/F/Z data blocks into storage. Band index advances on each "# lambda".
      integer, intent(in)  :: u
      integer, intent(out) :: status
      character(len=512) :: line, s, content
      integer :: ios, iband, i, is, ip, ii, c, irow, jcol
      real(wp) :: th, krow(4), frow(12), a1, a2, a3, zrow(16)

      status = 0;  iband = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         s = adjustl(line)
         if (len_trim(s) == 0) cycle
         if (s(1:1) /= '#') cycle
         content = adjustl(s(2:))

         if (starts_with(content, 'alignment profile:')) then
            call parse_profile(content)

         else if (starts_with(content, 'lambda')) then
            iband = iband + 1
            if (iband > scm_nband) then;  status = 2;  return;  end if
            scm_lambda(iband) = real_after_first_eq(content)

         else if (starts_with(content, 'Cext_tot')) then
            scm_cext_tot(iband) = real_after_first_eq(content)
         else if (starts_with(content, 'Csca_tot')) then
            scm_csca_tot(iband) = real_after_first_eq(content)
         else if (starts_with(content, 'Cext_ref')) then
            scm_cext_ref(iband) = real_after_first_eq(content)
         else if (starts_with(content, 'Csca_ref')) then
            scm_csca_ref(iband) = real_after_first_eq(content)

         else if (starts_with(content, 'K block')) then
            do i = 1, scm_nti
               read(u,*,iostat=ios) th, krow(1:4)
               if (ios /= 0) then;  status = 3;  return;  end if
               if (iband == 1) then
                  scm_theta_i(i) = th
               else if (abs(th - scm_theta_i(i)) > 1.0e-6_wp) then
                  status = 4;  return
               end if
               scm_cext_al(i, iband) = krow(1)
               scm_cpol_al(i, iband) = krow(2)
               scm_cbir_al(i, iband) = krow(3)
               scm_csca_al(i, iband) = krow(4)
            end do

         else if (starts_with(content, 'F block')) then
            do i = 1, scm_ntheta
               read(u,*,iostat=ios) th, frow(1:12)
               if (ios /= 0) then;  status = 3;  return;  end if
               if (iband == 1) then
                  scm_theta_ran(i) = th
               else if (abs(th - scm_theta_ran(i)) > 1.0e-6_wp) then
                  status = 4;  return
               end if
               scm_F_tot(i, 1:6, iband) = frow(1:6)
               scm_F_ref(i, 1:6, iband) = frow(7:12)
            end do

         else if (starts_with(content, 'Z block')) then
            do ii = 1, scm_nti
               do is = 1, scm_nts
                  do ip = 1, scm_nphi
                     read(u,*,iostat=ios) a1, a2, a3, zrow(1:16)
                     if (ios /= 0) then;  status = 3;  return;  end if
                     ! The Z block carries its own theta_s and phi grids (col 2
                     ! and col 3); build them from the first band.
                     if (iband == 1 .and. ii == 1 .and. ip == 1) scm_theta_s(is) = a2
                     if (iband == 1 .and. ii == 1 .and. is == 1) scm_phi(ip)     = a3
                     if (abs(a1 - scm_theta_i(ii)) > 1.0e-4_wp .or. &
                         abs(a2 - scm_theta_s(is)) > 1.0e-4_wp .or. &
                         abs(a3 - scm_phi(ip))     > 1.0e-4_wp) then
                        status = 4;  return
                     end if
                     do c = 1, 16
                        irow = (c-1)/4 + 1
                        jcol = mod(c-1, 4) + 1
                        scm_Z(ii, is, ip, irow, jcol, iband) = zrow(c)
                     end do
                  end do
               end do
            end do
         end if
      end do

      if (iband /= scm_nband) status = 2
   end subroutine read_bands


   subroutine allocate_storage()
      allocate(scm_lambda(scm_nband))
      allocate(scm_theta_i(scm_nti), scm_theta_s(scm_nts), scm_phi(scm_nphi))
      allocate(scm_theta_ran(scm_ntheta))
      allocate(scm_cos_theta_s(scm_nts))
      allocate(scm_cext_al(scm_nti, scm_nband), scm_cpol_al(scm_nti, scm_nband), &
               scm_cbir_al(scm_nti, scm_nband), scm_csca_al(scm_nti, scm_nband))
      allocate(scm_cext_tot(scm_nband), scm_csca_tot(scm_nband), &
               scm_cext_ref(scm_nband), scm_csca_ref(scm_nband))
      allocate(scm_F_tot(scm_ntheta, 6, scm_nband), scm_F_ref(scm_ntheta, 6, scm_nband))
      allocate(scm_Z(scm_nti, scm_nts, scm_nphi, 4, 4, scm_nband))
   end subroutine allocate_storage


   subroutine free_scatmat_aligned()
      if (allocated(scm_lambda))      deallocate(scm_lambda)
      if (allocated(scm_theta_i))     deallocate(scm_theta_i)
      if (allocated(scm_theta_s))     deallocate(scm_theta_s)
      if (allocated(scm_phi))         deallocate(scm_phi)
      if (allocated(scm_theta_ran))   deallocate(scm_theta_ran)
      if (allocated(scm_cos_theta_s)) deallocate(scm_cos_theta_s)
      if (allocated(scm_cext_al))     deallocate(scm_cext_al)
      if (allocated(scm_cpol_al))     deallocate(scm_cpol_al)
      if (allocated(scm_cbir_al))     deallocate(scm_cbir_al)
      if (allocated(scm_csca_al))     deallocate(scm_csca_al)
      if (allocated(scm_cext_tot))    deallocate(scm_cext_tot)
      if (allocated(scm_csca_tot))    deallocate(scm_csca_tot)
      if (allocated(scm_cext_ref))    deallocate(scm_cext_ref)
      if (allocated(scm_csca_ref))    deallocate(scm_csca_ref)
      if (allocated(scm_F_tot))       deallocate(scm_F_tot)
      if (allocated(scm_F_ref))       deallocate(scm_F_ref)
      if (allocated(scm_Z))           deallocate(scm_Z)
      scm_loaded = .false.
      scm_nband = 0;  scm_nti = 0;  scm_nts = 0;  scm_nphi = 0;  scm_ntheta = 0
      scm_profile_name = '';  scm_fmax = 0.0_wp
      scm_a_align = 0.0_wp;    scm_alpha = 0.0_wp
      scm_profile_mismatch = .false.
      scm_bytes = 0_int64
   end subroutine free_scatmat_aligned


   integer(int64) function storage_bytes() result(nb)
      integer(int64), parameter :: R8 = 8_int64
      nb = 0_int64
      nb = nb + int(scm_nband, int64) * R8                                    ! lambda
      nb = nb + int(scm_nti + scm_nts + scm_nphi + scm_ntheta + scm_nts, int64) * R8  ! grids + cos
      nb = nb + int(4*scm_nti, int64) * int(scm_nband, int64) * R8            ! K block
      nb = nb + 4_int64 * int(scm_nband, int64) * R8                          ! band scalars
      nb = nb + 2_int64 * int(scm_ntheta, int64) * 6_int64 * int(scm_nband, int64) * R8  ! F
      nb = nb + int(scm_nti, int64) * int(scm_nts, int64) * int(scm_nphi, int64) &
                * 16_int64 * int(scm_nband, int64) * R8                       ! Z
   end function storage_bytes


   ! ====================================================================
   ! Alignment-consistency guard
   ! ====================================================================

   subroutine alignment_matches_scatmat(f_max, a_align, alpha_align, tabulated, matched)
      ! Compare a requested alignment profile against the one the loaded aligned
      ! scattering table was integrated under. When no table is loaded there is
      ! nothing to conflict with, so matched = .true. A tabulated profile never
      ! matches the recorded analytic profile. Sets scm_profile_mismatch.
      !
      ! Rationale. The K and Z arrays were integrated over size once, under the
      ! recorded f_align profile; the sanctioned runtime variation is the
      ! scalar eta, not a change of profile. A different profile requires
      ! regenerating the table with run_scatmat_aligned.x profile=FILE.
      real(wp), intent(in)  :: f_max, a_align, alpha_align
      logical,  intent(in)  :: tabulated
      logical,  intent(out) :: matched

      if (.not. scm_loaded) then
         matched = .true.
         scm_profile_mismatch = .false.
         return
      end if

      if (tabulated) then
         matched = .false.
      else
         matched = rel_close(f_max,       scm_fmax)   .and. &
                   rel_close(a_align,     scm_a_align) .and. &
                   rel_close(alpha_align, scm_alpha)
      end if
      scm_profile_mismatch = .not. matched
   end subroutine alignment_matches_scatmat


   ! ====================================================================
   ! Query layer -- pure reads, thread-safe (no module state is written)
   ! ====================================================================

   subroutine scatmat_band(lambda_um, iband, exact)
      ! Nearest stored band to lambda_um, and whether it matched to BAND_TOL.
      ! Hot-path queries then take iband so no wavelength search runs per event.
      real(wp), intent(in)  :: lambda_um
      integer,  intent(out) :: iband
      logical,  intent(out) :: exact
      integer  :: k
      real(wp) :: d, dbest

      iband = 1;  dbest = abs(lambda_um - scm_lambda(1))
      do k = 2, scm_nband
         d = abs(lambda_um - scm_lambda(k))
         if (d < dbest) then
            dbest = d;  iband = k
         end if
      end do
      exact = (dbest <= BAND_TOL * abs(scm_lambda(iband)))
   end subroutine scatmat_band


   subroutine extinction_matrix_aligned(iband, theta_i, eta, kmat)
      ! kmat(4,4) [um^2/H]: the eta-scaled aligned extinction (attenuation)
      ! matrix at incidence theta_i, by linear interpolation in theta_i.
      !
      ! Standard aligned-spheroid form in the Mishchenko meridional basis
      ! (Q = Iv - Ih):
      !   diagonal   = eta * Cext_al
      !   K(1,2)=K(2,1) = eta * Cpol_al      (dichroism, IQ block)
      !   K(3,4)     = eta * Cbir_al         (birefringence, UV block)
      !   K(4,3)     = -eta * Cbir_al
      ! Sign convention: Cbir_al = 0.5 (Cre3 - Cre2) as recorded in the table
      ! header. The +Cbir in K(3,4) with -Cbir in K(4,3) is the antisymmetric
      ! circular-retardance block that rotates U into V on propagation; the RT
      ! rotates into the frame where Q is set by the projected field direction.
      !
      ! theta_i in [0,180] is folded to [0,90]: K is even under
      ! theta_i -> 180-theta_i because the forward propagation direction maps to
      ! itself under the equatorial reflection z -> -z of the oblate spheroid.
      integer,  intent(in)  :: iband
      real(wp), intent(in)  :: theta_i, eta
      real(wp), intent(out) :: kmat(4,4)
      integer  :: il
      real(wp) :: ti, t, ce, cp, cb

      ti = theta_i
      if (ti > 90.0_wp) ti = 180.0_wp - ti
      call bracket(scm_theta_i, scm_nti, ti, il, t)
      ce = (1.0_wp-t)*scm_cext_al(il, iband) + t*scm_cext_al(il+1, iband)
      cp = (1.0_wp-t)*scm_cpol_al(il, iband) + t*scm_cpol_al(il+1, iband)
      cb = (1.0_wp-t)*scm_cbir_al(il, iband) + t*scm_cbir_al(il+1, iband)

      kmat = 0.0_wp
      kmat(1,1) = eta*ce;  kmat(2,2) = eta*ce
      kmat(3,3) = eta*ce;  kmat(4,4) = eta*ce
      kmat(1,2) = eta*cp;  kmat(2,1) = eta*cp
      kmat(3,4) = eta*cb;  kmat(4,3) = -eta*cb
   end subroutine extinction_matrix_aligned


   subroutine mueller_matrix_aligned(iband, theta_i, theta_s, phi, z)
      ! z(4,4) [um^2 sr^-1 per H] at the REFERENCE alignment (eta=1); the caller
      ! scales by eta. Trilinear interpolation on (theta_i, theta_s, phi) with
      ! the header's symmetry reconstruction of the unstored ranges:
      !
      !   phi -> 360-phi        : off-diagonal 2x2 blocks flip sign.
      !   theta_i -> 180-theta_i: maps to theta_s -> 180-theta_s (phi unchanged)
      !                           with the same off-diagonal-block sign flip.
      !
      ! ORDER: interpolate first on the stored (folded) values, then apply the
      ! sign. The two reflections each contribute a constant +/-1 that is the
      ! SAME for every corner of the interpolation cell (it depends only on the
      ! query octant, not on position), so it factors out of the linear
      ! interpolation and the result is exact at grid nodes either way; applying
      ! it once to the interpolated matrix avoids touching eight corner blocks.
      integer,  intent(in)  :: iband
      real(wp), intent(in)  :: theta_i, theta_s, phi
      real(wp), intent(out) :: z(4,4)
      integer  :: ii, is, ip, i, j
      real(wp) :: ti, ts, ph, wi, ws, wp_, sgn
      real(wp) :: c000, c100, c010, c110, c001, c101, c011, c111

      ti = theta_i;  ts = theta_s;  ph = phi
      sgn = 1.0_wp

      ! Fold theta_i into [0,90] (with the coupled theta_s reflection).
      if (ti > 90.0_wp) then
         ti = 180.0_wp - ti
         ts = 180.0_wp - ts
         sgn = -sgn
      end if
      ! Fold phi into [0,180] via the phi mirror.
      ph = modulo(ph, 360.0_wp)
      if (ph > 180.0_wp) then
         ph = 360.0_wp - ph
         sgn = -sgn
      end if

      call bracket(scm_theta_i, scm_nti,  ti, ii, wi)
      call bracket(scm_theta_s, scm_nts,  ts, is, ws)
      call bracket(scm_phi,     scm_nphi, ph, ip, wp_)

      do j = 1, 4
         do i = 1, 4
            c000 = scm_Z(ii,   is,   ip,   i, j, iband)
            c100 = scm_Z(ii+1, is,   ip,   i, j, iband)
            c010 = scm_Z(ii,   is+1, ip,   i, j, iband)
            c110 = scm_Z(ii+1, is+1, ip,   i, j, iband)
            c001 = scm_Z(ii,   is,   ip+1, i, j, iband)
            c101 = scm_Z(ii+1, is,   ip+1, i, j, iband)
            c011 = scm_Z(ii,   is+1, ip+1, i, j, iband)
            c111 = scm_Z(ii+1, is+1, ip+1, i, j, iband)
            z(i,j) = (1.0_wp-wp_) * ( (1.0_wp-ws)*((1.0_wp-wi)*c000 + wi*c100) &
                                     +        ws *((1.0_wp-wi)*c010 + wi*c110) ) &
                   +        wp_  * ( (1.0_wp-ws)*((1.0_wp-wi)*c001 + wi*c101) &
                                     +        ws *((1.0_wp-wi)*c011 + wi*c111) )
         end do
      end do

      if (sgn < 0.0_wp) call flip_offdiagonal_blocks(z)
   end subroutine mueller_matrix_aligned


   subroutine mueller_matrix_random(iband, big_theta, f_tot, f_ref)
      ! The two six-element (11, 22, 33, 44, 12, 34) random-orientation matrices
      ! at scattering angle big_theta [deg], linear interpolation on the
      ! 1-degree Theta grid. Values are alpha1-normalized as stored; restore
      ! absolute um^2 sr^-1 per H with f_tot*scm_csca_tot, f_ref*scm_csca_ref.
      integer,  intent(in)  :: iband
      real(wp), intent(in)  :: big_theta
      real(wp), intent(out) :: f_tot(6), f_ref(6)
      integer  :: il, c
      real(wp) :: t

      call bracket(scm_theta_ran, scm_ntheta, big_theta, il, t)
      do c = 1, 6
         f_tot(c) = (1.0_wp-t)*scm_F_tot(il, c, iband) + t*scm_F_tot(il+1, c, iband)
         f_ref(c) = (1.0_wp-t)*scm_F_ref(il, c, iband) + t*scm_F_ref(il+1, c, iband)
      end do
   end subroutine mueller_matrix_random


   subroutine scattering_cross_sections(iband, theta_i, eta, csca_aligned, csca_unaligned)
      ! csca_aligned   = eta * Csca_al(theta_i)   [um^2/H], the fixed-orientation
      !                  aligned scattering cross section at incidence theta_i
      !                  (grid closure INT Z11 dOmega), folded to [0,90].
      ! csca_unaligned = Csca_tot - eta * Csca_ref [um^2/H], the random remainder
      !                  (theta_i-independent).
      integer,  intent(in)  :: iband
      real(wp), intent(in)  :: theta_i, eta
      real(wp), intent(out) :: csca_aligned, csca_unaligned
      integer  :: il
      real(wp) :: ti, t

      ti = theta_i
      if (ti > 90.0_wp) ti = 180.0_wp - ti
      call bracket(scm_theta_i, scm_nti, ti, il, t)
      csca_aligned   = eta * ((1.0_wp-t)*scm_csca_al(il, iband) + t*scm_csca_al(il+1, iband))
      csca_unaligned = scm_csca_tot(iband) - eta*scm_csca_ref(iband)
   end subroutine scattering_cross_sections


   ! ====================================================================
   ! Internal
   ! ====================================================================

   pure subroutine bracket(grid, n, x, ilo, t)
      ! Locate x in the strictly increasing grid(1:n): return ilo in [1,n-1] and
      ! the linear weight t in [0,1] for the segment grid(ilo)..grid(ilo+1),
      ! clamped at both ends. At any node x = grid(k), t is exactly 0 or 1, so an
      ! interpolation returns the stored node value bit-for-bit.
      real(wp), intent(in)  :: grid(:)
      integer,  intent(in)  :: n
      real(wp), intent(in)  :: x
      integer,  intent(out) :: ilo
      real(wp), intent(out) :: t
      integer :: lo, hi, mid

      if (x <= grid(1)) then
         ilo = 1;  t = 0.0_wp;  return
      end if
      if (x >= grid(n)) then
         ilo = n-1;  t = 1.0_wp;  return
      end if
      lo = 1;  hi = n
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (grid(mid) <= x) then
            lo = mid
         else
            hi = mid
         end if
      end do
      ilo = lo
      t = (x - grid(ilo)) / (grid(ilo+1) - grid(ilo))
   end subroutine bracket


   pure subroutine flip_offdiagonal_blocks(z)
      ! Negate the two off-diagonal 2x2 blocks (elements 13,14,23,24,31,32,41,42)
      ! -- the sign pattern of both the phi mirror and the equatorial mapping.
      real(wp), intent(inout) :: z(4,4)
      z(1,3) = -z(1,3);  z(1,4) = -z(1,4)
      z(2,3) = -z(2,3);  z(2,4) = -z(2,4)
      z(3,1) = -z(3,1);  z(3,2) = -z(3,2)
      z(4,1) = -z(4,1);  z(4,2) = -z(4,2)
   end subroutine flip_offdiagonal_blocks


   pure logical function rel_close(a, b)
      real(wp), intent(in) :: a, b
      real(wp) :: scale
      scale = max(abs(a), abs(b), tiny(1.0_wp))
      rel_close = (abs(a - b) <= PROFILE_TOL * scale)
   end function rel_close


   pure logical function starts_with(s, pre)
      character(len=*), intent(in) :: s, pre
      integer :: lp
      lp = len_trim(pre)
      starts_with = .false.
      if (lp <= len(s)) starts_with = (s(1:lp) == pre(1:lp))
   end function starts_with


   real(wp) function real_after_first_eq(line) result(val)
      ! Read a real from the text following the first '=' in line.
      character(len=*), intent(in) :: line
      integer :: p, ios
      val = 0.0_wp
      p = index(line, '=')
      if (p == 0 .or. p >= len(line)) return
      read(line(p+1:), *, iostat=ios) val
      if (ios /= 0) val = 0.0_wp
   end function real_after_first_eq


   subroutine parse_profile(content)
      ! content (after the leading '#') = "alignment profile: NAME  f_max = ...,
      ! a_align = ... um, alpha = ... (...)". Extract NAME and the three scalars.
      character(len=*), intent(in) :: content
      integer :: p, ios
      p = index(content, ':')
      scm_profile_name = ''
      if (p > 0 .and. p < len(content)) then
         read(content(p+1:), *, iostat=ios) scm_profile_name
         if (ios /= 0) scm_profile_name = ''
      end if
      scm_fmax    = real_after_key(content, 'f_max')
      scm_a_align = real_after_key(content, 'a_align')
      scm_alpha   = real_after_key(content, 'alpha')
   end subroutine parse_profile


   real(wp) function real_after_key(line, key) result(val)
      ! Read a real from the text following the first '=' at or after key.
      character(len=*), intent(in) :: line, key
      integer :: k, e, ios
      val = 0.0_wp
      k = index(line, key)
      if (k == 0) return
      e = index(line(k:), '=')
      if (e == 0) return
      read(line(k+e:), *, iostat=ios) val
      if (ios /= 0) val = 0.0_wp
   end function real_after_key


   subroutine gunzip_to(gz_file, out_file, ok)
      character(len=*), intent(in)  :: gz_file, out_file
      logical,          intent(out) :: ok
      integer :: estat, cstat
      call execute_command_line('gzip -dc "'//trim(gz_file)//'" > "'// &
                                trim(out_file)//'"', exitstat=estat, cmdstat=cstat)
      ok = (cstat == 0 .and. estat == 0)
   end subroutine gunzip_to


   subroutine discard_scratch(gz, path)
      logical,          intent(in) :: gz
      character(len=*), intent(in) :: path
      integer :: u, ios
      if (.not. gz) return
      open(newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close(u, status='delete')
   end subroutine discard_scratch

end module scatmat_aligned_mod
