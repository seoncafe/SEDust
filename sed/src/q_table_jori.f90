module q_table_jori_mod
   ! Orientation-resolved optics for the Draine & Hensley (2021) astrodust
   ! spheroid, plus the polarized cross sections and the grain-alignment
   ! efficiency they feed.
   !
   ! Input file: data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz
   !
   !   12 header lines, then free-format values in the order
   !     ((Q(jw,jr,jori), jw=0,1128), jr=0,168), jori=1,3
   !   written once for Q_ext, once for Q_abs, once for Q_sca.
   !   On disk each record holds the 169 sizes of one (jw, jori) pair, so
   !   the stream is read as 3 quantities x 3 orientations x 1129 records.
   !
   !   jori=1: k || a          (a = spheroid symmetry axis)
   !   jori=2: k perp a, E || a
   !   jori=3: k perp a, E perp a
   !
   ! The wavelength and size axes are NOT parsed out of the header; they are
   ! read from the companion grid files data/dielectric/DH21_wave and
   ! data/dielectric/DH21_aeff, which list the same nodes the table was
   ! computed on.
   !
   ! Derived quantities, for a grain whose symmetry axis is perpendicular to
   ! the line of sight and perfectly aligned:
   !
   !   Q_pol = 0.5 * (Q(jori=3) - Q(jori=2))      polarization cross section
   !   Q_ran = (Q(1) + Q(2) + Q(3)) / 3           random-orientation average
   !
   ! and a cross section follows from C = Q * pi * a_eff^2 with a_eff in cm.
   !
   ! gzip handling: the table ships compressed and Fortran cannot read a
   ! deflate stream, so the reader shells out to `gzip -dc` once, writes a
   ! scratch copy next to the caller's working directory, reads it, and
   ! deletes it. This keeps a fresh clone working with no manual setup step
   ! and leaves neither an untracked ~12 MB sibling in data/dielectric nor a
   ! dependency on a zlib binding. A path that does not end in `.gz` is
   ! opened directly.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: load_q_table_jori, falign_hd23
   public :: nj_lam, nj_aeff, lam_j, aeff_j
   public :: qext_j, qabs_j, qsca_j
   public :: qpol_ext, qpol_abs, qran_ext, qran_abs, qran_sca
   public :: A_ALIGN, ALPHA_ALIGN, FMAX_ALIGN

   integer, parameter :: wp = real64

   integer, parameter :: NA_DEF = 169    ! DH21 size grid length
   integer, parameter :: NW_DEF = 1129   ! DH21 wavelength grid length
   integer, parameter :: NHEAD  = 12     ! header lines in the Q table
   integer, parameter :: NORI   = 3      ! orientations stored

   ! Alignment efficiency, Hensley & Draine (2023) Table 1:
   !   f_align(a) = f_max / (1 + (a_align/a)**alpha_align)
   real(wp), parameter :: A_ALIGN     = 0.0749_wp   ! [um]
   real(wp), parameter :: ALPHA_ALIGN = 1.80_wp
   real(wp), parameter :: FMAX_ALIGN  = 1.00_wp

   integer  :: nj_lam = 0, nj_aeff = 0
   real(wp), allocatable :: lam_j(:), aeff_j(:)                ! grid axes [um]
   real(wp), allocatable :: qext_j(:,:,:), qabs_j(:,:,:), qsca_j(:,:,:)  ! (NLAM, NA, 3)
   real(wp), allocatable :: qpol_ext(:,:), qpol_abs(:,:)       ! (NLAM, NA)
   real(wp), allocatable :: qran_ext(:,:), qran_abs(:,:), qran_sca(:,:)

contains

   ! Alignment efficiency for a grain of effective radius a [um].
   !
   ! This analytic form is the canonical one for SEDust. Note that
   ! data/release/size_distribution.dat also ships an f_align column; the two
   ! agree to a median ratio of 1.0000 but differ by up to 0.32% at the small
   ! end. The formula is preferred because it is the published fit, it is
   ! defined for any radius grid rather than the tabulated one, and it keeps
   ! the alignment model in one place. Do not mix the two.
   pure function falign_hd23(a) result(f)
      real(wp), intent(in) :: a
      real(wp)             :: f
      if (a <= 0.0_wp) then
         f = 0.0_wp
      else
         f = FMAX_ALIGN / (1.0_wp + (A_ALIGN / a)**ALPHA_ALIGN)
      end if
   end function falign_hd23


   subroutine load_q_table_jori(q_file, wave_file, aeff_file, ok)
      ! Reads the orientation-resolved DH21 Q table and its grid axes, then
      ! fills the polarized and random-orientation combinations.
      !
      ! Optional ok follows load_q_table's convention: absent -> print and
      ! stop on any error; present -> return .false. with the module left
      ! unloaded (arrays freed, nj_lam = nj_aeff = 0) so a host can recover.
      character(len=*),  intent(in)  :: q_file, wave_file, aeff_file
      logical, optional, intent(out) :: ok

      integer  :: u, ios, iq, jori, jw, ja, i
      logical  :: gz, sub_ok
      real(wp) :: xextra
      real(wp), allocatable :: row(:)
      character(len=512)    :: read_path, line

      if (present(ok)) ok = .true.

      ! free_state() zeroes the grid counters, so set them only afterwards.
      call free_state()
      nj_aeff = NA_DEF
      nj_lam  = NW_DEF

      allocate(lam_j(nj_lam), aeff_j(nj_aeff), row(nj_aeff))
      allocate(qext_j(nj_lam, nj_aeff, NORI), qabs_j(nj_lam, nj_aeff, NORI), &
               qsca_j(nj_lam, nj_aeff, NORI))
      allocate(qpol_ext(nj_lam, nj_aeff), qpol_abs(nj_lam, nj_aeff))
      allocate(qran_ext(nj_lam, nj_aeff), qran_abs(nj_lam, nj_aeff), &
               qran_sca(nj_lam, nj_aeff))

      ! ---- grid axes -------------------------------------------------
      call read_grid(wave_file, nj_lam, lam_j, sub_ok)
      if (.not. sub_ok) then
         call bail('cannot read wavelength grid '//trim(wave_file))
         return
      end if
      call read_grid(aeff_file, nj_aeff, aeff_j, sub_ok)
      if (.not. sub_ok) then
         call bail('cannot read size grid '//trim(aeff_file))
         return
      end if

      do jw = 2, nj_lam
         if (lam_j(jw) <= lam_j(jw-1)) then
            call bail('lam_j not strictly increasing')
            return
         end if
      end do
      do ja = 2, nj_aeff
         if (aeff_j(ja) <= aeff_j(ja-1)) then
            call bail('aeff_j not strictly increasing')
            return
         end if
      end do

      ! ---- decompress if needed --------------------------------------
      gz = .false.
      i  = len_trim(q_file)
      if (i > 3) gz = (q_file(i-2:i) == '.gz')

      if (gz) then
         read_path = 'q_jori_scratch.dat'
         call gunzip_to(q_file, trim(read_path), sub_ok)
         if (.not. sub_ok) then
            ! The redirection creates the target before gzip can fail, so
            ! remove the empty file rather than leave it in the caller's
            ! working directory.
            call cleanup_scratch(.true., trim(read_path))
            call bail('gzip -dc failed on '//trim(q_file))
            return
         end if
      else
         read_path = q_file
      end if

      ! ---- read the table --------------------------------------------
      open(newunit=u, file=trim(read_path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call cleanup_scratch(gz, trim(read_path))
         call bail('cannot open '//trim(read_path))
         return
      end if

      do i = 1, NHEAD
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then
            close(u);  call cleanup_scratch(gz, trim(read_path))
            call bail('unexpected EOF in header')
            return
         end if
      end do

      ! Stream order: quantity (ext, abs, sca) outermost, then orientation,
      ! then one record of nj_aeff values for each wavelength.
      do iq = 1, 3
         do jori = 1, NORI
            do jw = 1, nj_lam
               read(u,*,iostat=ios) row(1:nj_aeff)
               if (ios /= 0) then
                  close(u);  call cleanup_scratch(gz, trim(read_path))
                  call bail('read error in Q block')
                  return
               end if
               do ja = 1, nj_aeff
                  if (.not. ieee_is_finite(row(ja))) then
                     close(u);  call cleanup_scratch(gz, trim(read_path))
                     call bail('non-finite Q value')
                     return
                  end if
               end do
               select case (iq)
               case (1);  qext_j(jw, :, jori) = row(1:nj_aeff)
               case (2);  qabs_j(jw, :, jori) = row(1:nj_aeff)
               case (3);  qsca_j(jw, :, jori) = row(1:nj_aeff)
               end select
            end do
         end do
      end do

      ! Reject a file that carries more than the expected payload.
      read(u,*,iostat=ios) xextra
      if (ios == 0) then
         close(u);  call cleanup_scratch(gz, trim(read_path))
         call bail('file has more data than the declared grid')
         return
      end if
      close(u)
      call cleanup_scratch(gz, trim(read_path))

      ! ---- derived combinations --------------------------------------
      qpol_ext = 0.5_wp * (qext_j(:,:,3) - qext_j(:,:,2))
      qpol_abs = 0.5_wp * (qabs_j(:,:,3) - qabs_j(:,:,2))
      qran_ext = (qext_j(:,:,1) + qext_j(:,:,2) + qext_j(:,:,3)) / 3.0_wp
      qran_abs = (qabs_j(:,:,1) + qabs_j(:,:,2) + qabs_j(:,:,3)) / 3.0_wp
      qran_sca = (qsca_j(:,:,1) + qsca_j(:,:,2) + qsca_j(:,:,3)) / 3.0_wp

      deallocate(row)

   contains

      subroutine bail(msg)
         ! Report an error the way the caller asked for: status flag, or
         ! message plus stop when no flag was supplied.
         character(len=*), intent(in) :: msg
         if (present(ok)) then
            call free_state();  ok = .false.
         else
            write(error_unit,'(a,a)') 'load_q_table_jori: ', msg
            stop 1
         end if
      end subroutine bail

   end subroutine load_q_table_jori


   subroutine free_state()
      ! Drop everything and mark the table unloaded.
      if (allocated(lam_j))    deallocate(lam_j)
      if (allocated(aeff_j))   deallocate(aeff_j)
      if (allocated(qext_j))   deallocate(qext_j)
      if (allocated(qabs_j))   deallocate(qabs_j)
      if (allocated(qsca_j))   deallocate(qsca_j)
      if (allocated(qpol_ext)) deallocate(qpol_ext)
      if (allocated(qpol_abs)) deallocate(qpol_abs)
      if (allocated(qran_ext)) deallocate(qran_ext)
      if (allocated(qran_abs)) deallocate(qran_abs)
      if (allocated(qran_sca)) deallocate(qran_sca)
      nj_lam = 0;  nj_aeff = 0
   end subroutine free_state


   subroutine read_grid(filename, n, arr, ok)
      ! DH21_wave / DH21_aeff: two title lines, then n free-format values.
      character(len=*), intent(in)  :: filename
      integer,          intent(in)  :: n
      real(wp),         intent(out) :: arr(:)
      logical,          intent(out) :: ok
      integer :: u, ios, i
      character(len=512) :: line

      ok = .false.
      arr(1:n) = 0.0_wp
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do i = 1, 2
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then
            close(u);  return
         end if
      end do
      read(u,*,iostat=ios) arr(1:n)
      close(u)
      if (ios /= 0) return
      do i = 1, n
         if (.not. ieee_is_finite(arr(i))) return
      end do
      ok = .true.
   end subroutine read_grid


   subroutine gunzip_to(gz_file, out_file, ok)
      ! Expand gz_file to out_file with the system gzip.
      character(len=*), intent(in)  :: gz_file, out_file
      logical,          intent(out) :: ok
      integer :: estat, cstat

      call execute_command_line('gzip -dc "'//trim(gz_file)//'" > "'// &
                                trim(out_file)//'"', &
                                exitstat=estat, cmdstat=cstat)
      ok = (cstat == 0 .and. estat == 0)
   end subroutine gunzip_to


   subroutine cleanup_scratch(gz, path)
      ! Remove the expanded copy, if we made one.
      logical,          intent(in) :: gz
      character(len=*), intent(in) :: path
      integer :: u, ios
      if (.not. gz) return
      open(newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close(u, status='delete')
   end subroutine cleanup_scratch

end module q_table_jori_mod
