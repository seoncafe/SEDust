module q_component_mod
   ! Reads a stored (lambda, a_eff) cross-section table of ONE grain
   ! population -- the DL07 and Zubko counterpart of the astrodust Q table
   ! that q_table.f90 reads.  sed/calc_qtable.x writes them into
   ! data/qtable/, both wavelength sets of every model:
   !
   !   q_dl07_{sil,gra,pah_neu,pah_ion}[_euv].dat.gz
   !   q_zubko_{sil,gra,pah}[_euv].dat.gz
   !
   ! Layout, as calc_qtable.x writes it: '#'-commented header, two lines of
   ! which carry the grid sizes,
   !
   !   # a_eff stride: every 1 of <NA>
   !   # lambda stride: every 1 of <NLAM>
   !
   ! then NLAM*NA rows, lambda-major with a_eff running fastest:
   !
   !   lambda[um]  a_eff[um]  Q_ext  Q_abs  Q_sca  albedo  g      (7 columns)
   !   lambda[um]  a_eff[um]  Q_abs                                (3 columns)
   !
   ! The three-column form is a population that does not scatter in the model
   ! -- the PAHs, whose cross sections are an absorption prescription rather
   ! than a Mie solution.  Q_sca and g come back zero for it, which is what
   ! the size integral must use.
   !
   ! gzip: these tables ship as plain text and are read directly.  A path that
   ! DOES end in '.gz' is still accepted -- Fortran cannot read a deflate
   ! stream, so it is expanded once with the system gzip into a scratch copy
   ! under TMPDIR (else /tmp), read, and deleted.  The scratch name carries the
   ! process id so concurrent ranks never collide.
   use, intrinsic :: iso_c_binding, only: c_int
   use constants, only: wp
   implicit none
   private
   public :: load_q_component

   interface
      function c_getpid() bind(c, name='getpid') result(pid)
         import :: c_int
         integer(c_int) :: pid
      end function c_getpid
   end interface

contains

   subroutine load_q_component(path, nlam, na, lam, aeff, Qabs, Qsca, gsca, ok, rho)
      ! On success lam(nlam), aeff(na) are the table's own axes and
      ! Qabs/Qsca/gsca are (nlam, na).  ok = .false. leaves everything
      ! deallocated and says why on stdout, so a builder can report it through
      ! its status argument instead of stopping an RT host.
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: nlam, na
      real(wp), allocatable, intent(out) :: lam(:), aeff(:)
      real(wp), allocatable, intent(out) :: Qabs(:,:), Qsca(:,:), gsca(:,:)
      logical,               intent(out) :: ok
      ! Solid density of the material [g/cm^3] when the table states one;
      ! left untouched otherwise, so a caller that has it already may pass a
      ! variable it has set.
      real(wp), optional,    intent(inout) :: rho

      character(len=512) :: read_path, line
      logical  :: gz
      integer  :: u, ios, i, jw, ja, ncol
      real(wp) :: l, a, Qe, Qa, Qs, alb, g

      ok = .false.;  nlam = 0;  na = 0

      gz = .false.
      i  = len_trim(path)
      if (i > 3) gz = (path(i-2:i) == '.gz')
      if (gz) then
         read_path = unique_scratch_path()
         call execute_command_line('gzip -dc "'//trim(path)//'" > "'// &
                                   trim(read_path)//'"', exitstat=ios)
         if (ios /= 0) then
            call discard_scratch_copy(.true., trim(read_path))
            write(*,'(a,a)') ' load_q_component: gzip -dc failed on ', trim(path)
            return
         end if
      else
         read_path = path
      end if

      open(newunit=u, file=trim(read_path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call discard_scratch_copy(gz, trim(read_path))
         write(*,'(a,a)') ' load_q_component: cannot open ', trim(path)
         return
      end if

      ! ---- header: the two grid sizes, and how wide the rows are ---------
      ncol = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (line(1:1) /= '#') then
            backspace(u)
            exit
         end if
         if (index(line,'a_eff stride') > 0) &
            read(line(index(line,'every 1 of')+10:), *, iostat=ios) na
         if (index(line,'lambda stride') > 0) &
            read(line(index(line,'every 1 of')+10:), *, iostat=ios) nlam
         if (index(line,'density [g/cm^3]') > 0 .and. present(rho)) &
            read(line(index(line,':')+1:), *, iostat=ios) rho
         if (index(line,'Q_sca') > 0) ncol = 7
         if (index(line,'Q_abs') > 0 .and. ncol == 0) ncol = 3
      end do
      if (nlam <= 0 .or. na <= 0 .or. ncol == 0) then
         close(u);  call discard_scratch_copy(gz, trim(read_path))
         write(*,'(a,a)') ' load_q_component: no grid sizes in the header of ', trim(path)
         return
      end if

      allocate(lam(nlam), aeff(na))
      allocate(Qabs(nlam,na), Qsca(nlam,na), gsca(nlam,na))
      Qsca = 0.0_wp;  gsca = 0.0_wp

      do jw = 1, nlam
         do ja = 1, na
            if (ncol == 7) then
               read(u,*,iostat=ios) l, a, Qe, Qa, Qs, alb, g
            else
               read(u,*,iostat=ios) l, a, Qa
               Qs = 0.0_wp;  g = 0.0_wp
            end if
            if (ios /= 0) then
               close(u);  call discard_scratch_copy(gz, trim(read_path))
               deallocate(lam, aeff, Qabs, Qsca, gsca)
               write(*,'(a,a)') ' load_q_component: short or malformed table ', trim(path)
               nlam = 0;  na = 0
               return
            end if
            if (jw == 1) aeff(ja) = a
            Qabs(jw,ja) = Qa;  Qsca(jw,ja) = Qs;  gsca(jw,ja) = g
         end do
         lam(jw) = l
      end do
      close(u)
      call discard_scratch_copy(gz, trim(read_path))
      ok = .true.
   end subroutine load_q_component


   function unique_scratch_path() result(p)
      ! Collision-free decompression target: TMPDIR (else /tmp) plus the
      ! process id, so concurrent MPI ranks never share the scratch name and a
      ! read-only launch directory is never written.
      character(len=512) :: p
      character(len=512) :: tmpdir
      character(len=32)  :: pidstr
      integer :: stat
      call get_environment_variable('TMPDIR', tmpdir, status=stat)
      if (stat /= 0 .or. len_trim(tmpdir) == 0) tmpdir = '/tmp'
      write(pidstr,'(i0)') int(c_getpid())
      p = trim(tmpdir)//'/q_component_'//trim(pidstr)//'.dat'
   end function unique_scratch_path


   subroutine discard_scratch_copy(gz, path)
      ! Remove the expanded copy, if we made one.
      logical,          intent(in) :: gz
      character(len=*), intent(in) :: path
      integer :: u, ios
      if (.not. gz) return
      open(newunit=u, file=path, status='old', iostat=ios)
      if (ios == 0) close(u, status='delete')
   end subroutine discard_scratch_copy

end module q_component_mod
