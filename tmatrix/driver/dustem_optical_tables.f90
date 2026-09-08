module dustem_optical_tables
   ! Readers for the two DustEM text layouts this directory needs:
   !
   !   LAMBDA.DAT      '#' comment lines, the number of wavelengths, then that
   !                   many wavelengths in microns (free format, wrapped).
   !
   !   Q_<gtype>.DAT   '#' comment lines, nsize, nsize radii in microns
   !   G_<gtype>.DAT   (wrapped 10 to a line), then one or more blocks, each
   !                   titled by a '####' comment line and holding one row per
   !                   LAMBDA.DAT wavelength of nsize values.  Q_ carries two
   !                   blocks (Q_abs then Q_sca); G_ carries one (<cos theta>).
   !
   ! sed/src/dustem_io.f90 reads the same layouts for the SED solver.  It is
   ! not reused here because it takes its working precision and a name-length
   ! parameter from the SED module chain, which this directory does not build;
   ! the two must be kept consistent, and any change to the layout belongs in
   ! both.  Nothing here writes: the G-file writer lives with the program that
   ! decides what goes in it.

   use, intrinsic :: iso_fortran_env, only: real64, error_unit
   implicit none
   private
   public :: read_dustem_wavelength_grid, read_dustem_sized_table

   integer, parameter :: wp = real64

contains

   subroutine read_dustem_wavelength_grid(path, nwave, lambda_um)
      character(len=*), intent(in) :: path
      integer, intent(out) :: nwave
      real(wp), allocatable, intent(out) :: lambda_um(:)
      integer :: u, ios

      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) call die('cannot open '//trim(path))
      call skip_comments(u)
      read(u,*,iostat=ios) nwave
      if (ios /= 0 .or. nwave < 2) call die('bad wavelength count in '//trim(path))
      allocate(lambda_um(nwave))
      call skip_comments(u)
      read(u,*,iostat=ios) lambda_um
      if (ios /= 0) call die('cannot read the wavelengths of '//trim(path))
      close(u)
   end subroutine read_dustem_wavelength_grid


   subroutine read_dustem_sized_table(path, nwave, nsize, a_um, block, nblock)
      !! Reads the radii and up to nblock lambda-by-size blocks.  nblock is
      !! the number requested on entry and the number actually read on exit.
      character(len=*), intent(in) :: path
      integer, intent(in) :: nwave
      integer, intent(out) :: nsize
      real(wp), allocatable, intent(out) :: a_um(:)
      real(wp), allocatable, intent(out) :: block(:,:,:)   ! (nwave, nsize, nblock)
      integer, intent(inout) :: nblock
      integer :: u, ios, ib, k, nb_want

      nb_want = nblock
      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) call die('cannot open '//trim(path))
      call skip_comments(u)
      read(u,*,iostat=ios) nsize
      if (ios /= 0 .or. nsize < 2) call die('bad size count in '//trim(path))
      allocate(a_um(nsize))
      call skip_comments(u)
      read(u,*,iostat=ios) a_um
      if (ios /= 0) call die('cannot read the radii of '//trim(path))

      allocate(block(nwave, nsize, nb_want))
      nblock = 0
      do ib = 1, nb_want
         call skip_comments(u)
         do k = 1, nwave
            read(u,*,iostat=ios) block(k, 1:nsize, ib)
            if (ios /= 0) call die('short block in '//trim(path))
         end do
         nblock = ib
      end do
      close(u)
   end subroutine read_dustem_sized_table


   subroutine skip_comments(u)
      !! Leaves the file positioned on the next line that is not blank and
      !! does not begin with '#'.
      integer, intent(in) :: u
      character(len=512) :: line
      integer :: ios
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) return
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         backspace(u)
         return
      end do
   end subroutine skip_comments


   subroutine die(msg)
      character(len=*), intent(in) :: msg
      write(error_unit,'(a,a)') ' dustem_optical_tables: ', msg
      stop 1
   end subroutine die

end module dustem_optical_tables
