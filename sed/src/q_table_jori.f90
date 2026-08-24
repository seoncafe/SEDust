module q_table_jori_mod
   ! Orientation-resolved optics for the Draine & Hensley (2021) astrodust
   ! spheroid, plus the polarized cross sections and the grain-alignment
   ! efficiency they feed.
   !
   ! Input file: data/astrodust/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz
   !
   !   12 header lines, then free-format values in the order
   !     ((Q(jw,jr,jori), jw=0,1128), jr=0,168), jori=1,3
   !   written once for Q_ext, once for Q_abs, once for Q_sca, and
   !   optionally a 4th time for Q_re (the real-part forward-amplitude twin
   !   of Q_ext).  On disk each record holds the 169 sizes of one (jw, jori)
   !   pair, so the stream is read as 3 (or 4) quantities x 3 orientations
   !   x 1129 records.  A 3-block table without Q_re still loads: has_bir is
   !   then .false. and qbir_ext is left unallocated.
   !
   !   jori=1: k || a          (a = spheroid symmetry axis)
   !   jori=2: k perp a, E || a
   !   jori=3: k perp a, E perp a
   !
   ! The wavelength and size axes are NOT parsed out of the header; they are
   ! read from the companion grid files data/astrodust/DH21_wave and
   ! data/astrodust/DH21_aeff, which list the same nodes the table was
   ! computed on.
   !
   ! Derived quantities, for a grain whose symmetry axis is perpendicular to
   ! the line of sight and perfectly aligned:
   !
   !   Q_pol = 0.5 * (Q(jori=3) - Q(jori=2))      polarization cross section
   !   Q_ran = (Q(1) + Q(2) + Q(3)) / 3           random-orientation average
   !
   ! and, when the optional Q_re block is present (has_bir = .true.),
   !
   !   Q_bir = 0.5 * (Q_re(jori=3) - Q_re(jori=2))  birefringence (U<->V phase
   !                                                retardation on propagation)
   !
   ! and a cross section follows from C = Q * pi * a_eff^2 with a_eff in cm.
   !
   ! THE EXTREME-ULTRAVIOLET COMPANION TABLE.  load_q_table_jori_euv reads a
   ! second file of exactly the same stream format computed on the wavelength
   ! axis data/astrodust/DH21_wave_euv, which runs from 0.0124 um (100 eV) up
   ! to the 0.0912 um (13.6 eV) node where the table above starts.  It is a
   ! separate entry point with its own module state (nj_*_euv, lam_j_euv,
   ! qpol_*_euv, qbir_ext_euv) because the two tables have different grid
   ! lengths and are loaded independently; the reader itself is shared
   ! (read_jori_stream) so the two go through identical parsing and validation.
   ! Its grid lengths come from the axis files rather than a compiled-in
   ! constant.  Nothing about load_q_table_jori changes when it is used.
   !
   ! gzip handling: the table ships compressed and Fortran cannot read a
   ! deflate stream, so the reader shells out to `gzip -dc` once, writes a
   ! scratch copy into TMPDIR (or /tmp), reads it, and deletes it. The scratch
   ! name carries the process id (q_jori_<pid>.dat), so concurrent MPI ranks
   ! never share it and a read-only launch directory is never written. This
   ! keeps a fresh clone working with no manual setup step and avoids a
   ! dependency on a zlib binding. A path that does not end in `.gz` is
   ! opened directly.
   use size_dist_mod, only: falign_hd23, falign_powerlaw, &
                            A_ALIGN, ALPHA_ALIGN, FMAX_ALIGN

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_c_binding,   only: c_int
   implicit none
   private

   interface
      ! POSIX getpid, used only to name a collision-free decompression scratch
      ! file so concurrent MPI ranks never write the same path. Local copy of
      ! the same helper in scatmat_aligned_mod; this repo keeps sibling copies.
      function c_getpid() bind(c, name="getpid") result(pid)
         import :: c_int
         integer(c_int) :: pid
      end function c_getpid
   end interface

   ! The HD23 alignment law and its constants live in size_dist_mod with the
   ! rest of the HD23 size-distribution fit; re-exported here for the
   ! polarized-optics callers.
   public :: load_q_table_jori, falign_hd23, falign_powerlaw
   ! The raw orientation-resolved arrays, before the dichroic and
   ! random-orientation combinations are formed from them.  Exposed so that
   ! sed/calc_polarized_optics.x can store what the table actually holds rather
   ! than a quantity derived from it.
   public :: read_jori_stream
   public :: nj_lam, nj_aeff, lam_j, aeff_j
   public :: qpol_ext, qpol_abs, qran_ext, qran_abs, qran_sca
   public :: qbir_ext, has_bir
   public :: A_ALIGN, ALPHA_ALIGN, FMAX_ALIGN
   ! Extreme-ultraviolet companion table (see the note above).
   public :: load_q_table_jori_euv, free_q_table_jori_euv
   public :: nj_lam_euv, nj_aeff_euv, lam_j_euv, aeff_j_euv
   public :: qpol_ext_euv, qpol_abs_euv, qbir_ext_euv, has_bir_euv

   integer, parameter :: wp = real64

   integer, parameter :: NA_DEF = 169    ! DH21 size grid length
   ! The wavelength grid length is NOT written down here.  run_q_jori.f90 takes
   ! its axis from a grid file and its length with it, so a table's length is a
   ! property of the axis it was swept on rather than of this module, and the
   ! orientation-resolved table need not match the scalar one: the shipped
   ! polarized table is the 1129-point DH21 axis, 0.0912-39810 um, while the
   ! scalar Q table now runs to 1.0e-4 um on 1762 points.  The length is
   ! counted out of the axis file the caller names, exactly as the EUV
   ! companion reader below already does.
   integer, parameter :: NHEAD  = 12     ! header lines in the Q table
   integer, parameter :: NORI   = 3      ! orientations stored


   integer  :: nj_lam = 0, nj_aeff = 0
   real(wp), allocatable :: lam_j(:), aeff_j(:)                ! grid axes [um]
   real(wp), allocatable :: qpol_ext(:,:), qpol_abs(:,:)       ! (NLAM, NA)
   real(wp), allocatable :: qran_ext(:,:), qran_abs(:,:), qran_sca(:,:)
   ! Birefringence 0.5*(Q_re(jori=3)-Q_re(jori=2)), formed from the optional
   ! 4th (Q_re) block when present.  has_bir is .false. for an older 3-block
   ! table, in which case qbir_ext is left unallocated.
   real(wp), allocatable :: qbir_ext(:,:)                      ! (NLAM, NA)
   logical  :: has_bir = .false.

   ! Extreme-ultraviolet companion table, 0.0124-0.0912 um.  Same quantities,
   ! its own grid, loaded and freed independently of the arrays above.
   integer  :: nj_lam_euv = 0, nj_aeff_euv = 0
   real(wp), allocatable :: lam_j_euv(:), aeff_j_euv(:)             ! [um]
   real(wp), allocatable :: qpol_ext_euv(:,:), qpol_abs_euv(:,:)    ! (NLAM_euv, NA)
   real(wp), allocatable :: qbir_ext_euv(:,:)                       ! (NLAM_euv, NA)
   logical  :: has_bir_euv = .false.

contains



   subroutine load_q_table_jori(q_file, wave_file, aeff_file, ok)
      ! Reads the orientation-resolved DH21 Q table and its grid axes, then
      ! fills the polarized and random-orientation combinations.
      !
      ! Optional ok follows load_q_table's convention: absent -> print and
      ! stop on any error; present -> return .false. with the module left
      ! unloaded (arrays freed, nj_lam = nj_aeff = 0) so a host can recover.
      character(len=*),  intent(in)  :: q_file, wave_file, aeff_file
      logical, optional, intent(out) :: ok

      logical  :: sub_ok, bir
      integer  :: nw_file
      real(wp), allocatable :: qext_j(:,:,:), qabs_j(:,:,:), qsca_j(:,:,:)  ! (NLAM, NA, 3)
      real(wp), allocatable :: qre_j(:,:,:)      ! (NLAM, NA, 3), optional 4th block
      character(len=512)    :: msg

      if (present(ok)) ok = .true.

      ! free_state() zeroes the grid counters, so set them only afterwards.
      call free_state()

      nw_file = count_grid_values(wave_file)
      if (nw_file < 2) then
         if (present(ok)) then
            call free_state();  ok = .false.
         else
            write(error_unit,'(a,a)') &
               'load_q_table_jori: cannot count the wavelength axis ', trim(wave_file)
            stop 1
         end if
         return
      end if

      call read_jori_stream(q_file, wave_file, aeff_file, NA_DEF, nw_file, &
                            lam_j, aeff_j, qext_j, qabs_j, qsca_j, qre_j, &
                            bir, sub_ok, msg)
      if (.not. sub_ok) then
         if (present(ok)) then
            call free_state();  ok = .false.
         else
            write(error_unit,'(a,a)') 'load_q_table_jori: ', trim(msg)
            stop 1
         end if
         return
      end if

      nj_aeff = NA_DEF
      nj_lam  = nw_file

      ! ---- derived combinations --------------------------------------
      allocate(qpol_ext(nj_lam, nj_aeff), qpol_abs(nj_lam, nj_aeff))
      allocate(qran_ext(nj_lam, nj_aeff), qran_abs(nj_lam, nj_aeff), &
               qran_sca(nj_lam, nj_aeff))
      qpol_ext = 0.5_wp * (qext_j(:,:,3) - qext_j(:,:,2))
      qpol_abs = 0.5_wp * (qabs_j(:,:,3) - qabs_j(:,:,2))
      qran_ext = (qext_j(:,:,1) + qext_j(:,:,2) + qext_j(:,:,3)) / 3.0_wp
      qran_abs = (qabs_j(:,:,1) + qabs_j(:,:,2) + qabs_j(:,:,3)) / 3.0_wp
      qran_sca = (qsca_j(:,:,1) + qsca_j(:,:,2) + qsca_j(:,:,3)) / 3.0_wp

      ! Birefringence, only from a 4-block table.
      has_bir = bir
      if (bir) then
         allocate(qbir_ext(nj_lam, nj_aeff))
         qbir_ext = 0.5_wp * (qre_j(:,:,3) - qre_j(:,:,2))
      end if

      ! The three orientations enter the optics only through Q_pol, Q_ran and
      ! the birefringence, so once those are formed the orientation-resolved
      ! arrays are spent; being locals, their ~13.7 MB go back on return rather
      ! than staying resident until unload.
   end subroutine load_q_table_jori


   subroutine load_q_table_jori_euv(q_file, wave_file, aeff_file, ok)
      ! Reads the extreme-ultraviolet companion table (0.0124-0.0912 um) and
      ! fills the polarized combinations on its own grid.  The wavelength and
      ! size axes are counted out of the axis files, so the EUV grid length is
      ! a property of the data rather than a compiled-in constant.
      !
      ! Only the polarized channels are formed: the random-orientation optics
      ! of this band come from the volume-equivalent Mie sphere in the SED
      ! model (q_astrodust_mod), not from here.
      !
      ! ok follows load_q_table_jori: absent -> message + stop; present ->
      ! .false. with the EUV state left unloaded.
      character(len=*),  intent(in)  :: q_file, wave_file, aeff_file
      logical, optional, intent(out) :: ok

      integer  :: na, nw
      logical  :: sub_ok, bir
      real(wp), allocatable :: qext_j(:,:,:), qabs_j(:,:,:), qsca_j(:,:,:)
      real(wp), allocatable :: qre_j(:,:,:)
      character(len=512)    :: msg

      if (present(ok)) ok = .true.
      call free_q_table_jori_euv()

      nw = count_grid_values(wave_file)
      na = count_grid_values(aeff_file)
      if (nw < 2 .or. na < 2) then
         call bail_euv('cannot count the EUV grid axes '//trim(wave_file)// &
                       ' / '//trim(aeff_file))
         return
      end if

      call read_jori_stream(q_file, wave_file, aeff_file, na, nw, &
                            lam_j_euv, aeff_j_euv, qext_j, qabs_j, qsca_j, &
                            qre_j, bir, sub_ok, msg)
      if (.not. sub_ok) then
         call bail_euv(trim(msg))
         return
      end if

      nj_lam_euv  = nw
      nj_aeff_euv = na
      allocate(qpol_ext_euv(nw, na), qpol_abs_euv(nw, na))
      qpol_ext_euv = 0.5_wp * (qext_j(:,:,3) - qext_j(:,:,2))
      qpol_abs_euv = 0.5_wp * (qabs_j(:,:,3) - qabs_j(:,:,2))
      has_bir_euv = bir
      if (bir) then
         allocate(qbir_ext_euv(nw, na))
         qbir_ext_euv = 0.5_wp * (qre_j(:,:,3) - qre_j(:,:,2))
      end if

   contains

      subroutine bail_euv(m)
         character(len=*), intent(in) :: m
         if (present(ok)) then
            call free_q_table_jori_euv();  ok = .false.
         else
            write(error_unit,'(a,a)') 'load_q_table_jori_euv: ', m
            stop 1
         end if
      end subroutine bail_euv

   end subroutine load_q_table_jori_euv


   subroutine read_jori_stream(q_file, wave_file, aeff_file, na, nw, &
                               lam_out, aeff_out, qe, qa, qs, qr, bir, ok, msg)
      ! Shared parser for both orientation-resolved tables: reads the two grid
      ! axes, checks them for strict monotonicity, decompresses the Q file if
      ! it is gzipped, skips the NHEAD header lines and reads the three
      ! mandatory (Q_ext, Q_abs, Q_sca) blocks plus the optional 4th (Q_re)
      ! block in (iq, jori, jw) stream order.
      !
      ! bir = .true. iff the 4th block was present, in which case qr is
      ! allocated; otherwise qr is left unallocated.  ok = .false. leaves every
      ! output array unallocated and puts the reason in msg.
      character(len=*), intent(in)  :: q_file, wave_file, aeff_file
      integer,          intent(in)  :: na, nw
      real(wp), allocatable, intent(out) :: lam_out(:), aeff_out(:)
      real(wp), allocatable, intent(out) :: qe(:,:,:), qa(:,:,:), qs(:,:,:), qr(:,:,:)
      logical,          intent(out) :: bir, ok
      character(len=*), intent(out) :: msg

      integer  :: u, ios, iq, jori, jw, ja, i
      logical  :: gz, sub_ok
      real(wp) :: xextra
      real(wp), allocatable :: row(:)
      character(len=512)    :: read_path, line

      ok  = .true.
      bir = .false.
      msg = ''

      allocate(lam_out(nw), aeff_out(na), row(na))
      allocate(qe(nw, na, NORI), qa(nw, na, NORI), qs(nw, na, NORI))

      ! ---- grid axes -------------------------------------------------
      call read_grid(wave_file, nw, lam_out, sub_ok)
      if (.not. sub_ok) then
         call fail('cannot read wavelength grid '//trim(wave_file))
         return
      end if
      call read_grid(aeff_file, na, aeff_out, sub_ok)
      if (.not. sub_ok) then
         call fail('cannot read size grid '//trim(aeff_file))
         return
      end if

      do jw = 2, nw
         if (lam_out(jw) <= lam_out(jw-1)) then
            call fail('lam_j not strictly increasing')
            return
         end if
      end do
      do ja = 2, na
         if (aeff_out(ja) <= aeff_out(ja-1)) then
            call fail('aeff_j not strictly increasing')
            return
         end if
      end do

      ! ---- decompress if needed --------------------------------------
      gz = .false.
      i  = len_trim(q_file)
      if (i > 3) gz = (q_file(i-2:i) == '.gz')

      if (gz) then
         read_path = unique_scratch_path()
         call gunzip_to(q_file, trim(read_path), sub_ok)
         if (.not. sub_ok) then
            ! The redirection creates the target before gzip can fail, so
            ! remove the empty file rather than leave it behind in TMPDIR.
            call discard_scratch_copy(.true., trim(read_path))
            call fail('gzip -dc failed on '//trim(q_file))
            return
         end if
      else
         read_path = q_file
      end if

      ! ---- read the table --------------------------------------------
      open(newunit=u, file=trim(read_path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call discard_scratch_copy(gz, trim(read_path))
         call fail('cannot open '//trim(read_path))
         return
      end if

      do i = 1, NHEAD
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then
            close(u);  call discard_scratch_copy(gz, trim(read_path))
            call fail('unexpected EOF in header')
            return
         end if
      end do

      ! Stream order: quantity (ext, abs, sca) outermost, then orientation,
      ! then one record of na values for each wavelength.
      do iq = 1, 3
         do jori = 1, NORI
            do jw = 1, nw
               read(u,*,iostat=ios) row(1:na)
               if (ios /= 0) then
                  close(u);  call discard_scratch_copy(gz, trim(read_path))
                  call fail('read error in Q block')
                  return
               end if
               do ja = 1, na
                  if (.not. ieee_is_finite(row(ja))) then
                     close(u);  call discard_scratch_copy(gz, trim(read_path))
                     call fail('non-finite Q value')
                     return
                  end if
               end do
               select case (iq)
               case (1);  qe(jw, :, jori) = row(1:na)
               case (2);  qa(jw, :, jori) = row(1:na)
               case (3);  qs(jw, :, jori) = row(1:na)
               end select
            end do
         end do
      end do

      ! ---- optional 4th block: Q_re (birefringence twin) -------------
      ! A 4-block table carries a further Q_re block in the same
      ! (jori, jw) nesting; its real-part forward amplitude gives the
      ! birefringence 0.5*(Qre3-Qre2).  An older 3-block table hits EOF right
      ! here: leave qr unallocated and bir = .false.
      read(u,*,iostat=ios) row(1:na)
      if (is_iostat_end(ios)) then
         bir = .false.
      else if (ios /= 0) then
         close(u);  call discard_scratch_copy(gz, trim(read_path))
         call fail('read error probing the Q_re block')
         return
      else
         allocate(qr(nw, na, NORI))
         ! The record just read is (jori = 1, jw = 1).
         if (.not. row_is_finite(row, na)) then
            close(u);  call discard_scratch_copy(gz, trim(read_path))
            call fail('non-finite Q_re value')
            return
         end if
         qr(1, :, 1) = row(1:na)
         do jori = 1, NORI
            do jw = 1, nw
               if (jori == 1 .and. jw == 1) cycle
               read(u,*,iostat=ios) row(1:na)
               if (ios /= 0) then
                  close(u);  call discard_scratch_copy(gz, trim(read_path))
                  call fail('read error in Q_re block')
                  return
               end if
               if (.not. row_is_finite(row, na)) then
                  close(u);  call discard_scratch_copy(gz, trim(read_path))
                  call fail('non-finite Q_re value')
                  return
               end if
               qr(jw, :, jori) = row(1:na)
            end do
         end do
         bir = .true.
      end if

      ! Reject a file that carries more than the expected payload, checked after
      ! whichever block was last (the 3rd for an old table, the 4th otherwise).
      read(u,*,iostat=ios) xextra
      if (ios == 0) then
         close(u);  call discard_scratch_copy(gz, trim(read_path))
         call fail('file has more data than the declared grid')
         return
      end if
      close(u)
      call discard_scratch_copy(gz, trim(read_path))

      deallocate(row)

   contains

      subroutine fail(m)
         ! Drop every partially filled output and report the reason.
         character(len=*), intent(in) :: m
         ok  = .false.
         bir = .false.
         msg = m
         if (allocated(lam_out))  deallocate(lam_out)
         if (allocated(aeff_out)) deallocate(aeff_out)
         if (allocated(qe)) deallocate(qe)
         if (allocated(qa)) deallocate(qa)
         if (allocated(qs)) deallocate(qs)
         if (allocated(qr)) deallocate(qr)
         if (allocated(row)) deallocate(row)
      end subroutine fail

      logical function row_is_finite(v, n)
         ! .true. iff v(1:n) are all finite.
         real(wp), intent(in) :: v(:)
         integer,  intent(in) :: n
         integer :: k
         row_is_finite = .true.
         do k = 1, n
            if (.not. ieee_is_finite(v(k))) then
               row_is_finite = .false.
               return
            end if
         end do
      end function row_is_finite

   end subroutine read_jori_stream


   integer function count_grid_values(filename) result(n)
      ! Number of whitespace-separated values in a grid-axis file after its two
      ! title lines.  Returns 0 when the file cannot be opened or holds none.
      !
      ! The record is consumed in fixed-size chunks with a non-advancing read
      ! rather than in one plain read: a grid file may put its whole axis on a
      ! single line (data/astrodust/DH21_wave puts all 1129 values on an
      ! 11289-character one), and a plain read into a fixed buffer would
      ! silently drop everything past the buffer and undercount. in_tok carries
      ! across chunks so a token straddling a boundary is counted once.
      character(len=*), intent(in) :: filename
      integer :: u, ios, i, L
      character(len=4096) :: line
      logical :: in_tok
      n = 0
      open(newunit=u, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      read(u,'(a)',iostat=ios) line
      read(u,'(a)',iostat=ios) line
      in_tok = .false.
      do
         read(u,'(a)',advance='no',size=L,iostat=ios) line
         do i = 1, L
            if (line(i:i) == ' ' .or. line(i:i) == char(9)) then
               in_tok = .false.
            else
               if (.not. in_tok) n = n + 1
               in_tok = .true.
            end if
         end do
         if (is_iostat_eor(ios)) then
            in_tok = .false.            ! end of record ends any open token
         else if (ios /= 0) then
            exit                        ! end of file
         end if
      end do
      close(u)
   end function count_grid_values


   subroutine free_q_table_jori_euv()
      ! Drop the EUV companion table and mark it unloaded.
      if (allocated(lam_j_euv))    deallocate(lam_j_euv)
      if (allocated(aeff_j_euv))   deallocate(aeff_j_euv)
      if (allocated(qpol_ext_euv)) deallocate(qpol_ext_euv)
      if (allocated(qpol_abs_euv)) deallocate(qpol_abs_euv)
      if (allocated(qbir_ext_euv)) deallocate(qbir_ext_euv)
      has_bir_euv = .false.
      nj_lam_euv = 0;  nj_aeff_euv = 0
   end subroutine free_q_table_jori_euv


   subroutine free_state()
      ! Drop everything and mark the table unloaded.
      if (allocated(lam_j))    deallocate(lam_j)
      if (allocated(aeff_j))   deallocate(aeff_j)
      if (allocated(qpol_ext)) deallocate(qpol_ext)
      if (allocated(qpol_abs)) deallocate(qpol_abs)
      if (allocated(qran_ext)) deallocate(qran_ext)
      if (allocated(qran_abs)) deallocate(qran_abs)
      if (allocated(qran_sca)) deallocate(qran_sca)
      if (allocated(qbir_ext)) deallocate(qbir_ext)
      has_bir = .false.
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


   function unique_scratch_path() result(p)
      ! Collision-free decompression target: TMPDIR (else /tmp) plus the process
      ! id, so concurrent MPI ranks never share the scratch name and a read-only
      ! launch directory is never written.
      character(len=512) :: p
      character(len=512) :: tmpdir
      character(len=32)  :: pidstr
      integer :: stat
      call get_environment_variable('TMPDIR', tmpdir, status=stat)
      if (stat /= 0 .or. len_trim(tmpdir) == 0) tmpdir = '/tmp'
      write(pidstr,'(i0)') int(c_getpid())
      p = trim(tmpdir)//'/q_jori_'//trim(pidstr)//'.dat'
   end function unique_scratch_path


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


   subroutine discard_scratch_copy(gz, path)
      ! Remove the expanded copy, if we made one.
      logical,          intent(in) :: gz
      character(len=*), intent(in) :: path
      integer :: u, ios
      if (.not. gz) return
      open(newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close(u, status='delete')
   end subroutine discard_scratch_copy

end module q_table_jori_mod
