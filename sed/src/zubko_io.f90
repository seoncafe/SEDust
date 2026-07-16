module zubko_io
   ! Readers for the Zubko, Dwek & Arendt (2004) BARE-GR-S dust model as
   ! distributed with the Camps et al. (2015) RT benchmark
   ! (SHG_Benchmark/DustModel). Three pieces, each a separate file format:
   !
   !   1. Size distribution -- ZDA log-polynomial FORMULA in an INI-style
   !      config (ZDA_BARE_GR_S_Config.dat), or the tabulated SzDist files.
   !   2. Optics -- DustEM/Zubko Q-tables (one block for each radius).
   !   3. Calorimetry -- specific enthalpy/heat-capacity tables.
   !
   ! This module currently implements the config parser + the ZDA size
   ! distribution formula (the user's preferred path); the optics and
   ! calorimetry readers and build_zubko follow.
   !
   ! ZDA size distribution (per component):
   !   log10 g(a) = c0 + b0*log10(a)
   !                   - b1*|log10(a/a1)|^m1 - b2*|log10(a/a2)|^m2
   !                   - b3*|a - a3|^m3       - b4*|a - a4|^m4
   !   f(a) = A * g(a),    a in micron,  [g] = um^-1
   ! Missing (NULL) terms have b_k = 0 and drop out.
   use constants, only: wp
   implicit none
   private
   public :: zda_comp_t, read_zda_config, zda_gofa, ZDA_MAXCOMP
   public :: read_zubko_optics, read_zubko_calor, read_dnda_table

   integer, parameter :: ZDA_MAXCOMP = 8

   type :: zda_comp_t
      character(len=32) :: label   = ''     ! e.g. 'PAH','Graphite','Silicate' (from comment)
      character(len=64) :: xsec    = ''      ! Cross Sections= file stem
      character(len=64) :: calor   = ''      ! Calorimetry= file
      real(wp) :: A     = 0.0_wp
      real(wp) :: a_min = 0.0_wp, a_max = 0.0_wp
      real(wp) :: c0 = 0.0_wp, b0 = 0.0_wp
      real(wp) :: b1 = 0.0_wp, a1 = 1.0_wp, m1 = 1.0_wp
      real(wp) :: b2 = 0.0_wp, a2 = 1.0_wp, m2 = 1.0_wp
      real(wp) :: b3 = 0.0_wp, a3 = 1.0_wp, m3 = 1.0_wp
      real(wp) :: b4 = 0.0_wp, a4 = 1.0_wp, m4 = 1.0_wp
   end type zda_comp_t

contains

   ! Evaluate f(a) = A*g(a) [um^-1 * A's units] for a component at a [um].
   ! Returns 0 outside [a_min, a_max].
   pure function zda_gofa(c, a_um) result(f)
      type(zda_comp_t), intent(in) :: c
      real(wp),         intent(in) :: a_um
      real(wp) :: f, lg, la
      if (a_um < c%a_min .or. a_um > c%a_max) then
         f = 0.0_wp; return
      end if
      la = log10(a_um)
      lg = c%c0 + c%b0*la
      if (c%b1 /= 0.0_wp) lg = lg - c%b1 * abs(log10(a_um/c%a1))**c%m1
      if (c%b2 /= 0.0_wp) lg = lg - c%b2 * abs(log10(a_um/c%a2))**c%m2
      if (c%b3 /= 0.0_wp) lg = lg - c%b3 * abs(a_um - c%a3)**c%m3
      if (c%b4 /= 0.0_wp) lg = lg - c%b4 * abs(a_um - c%a4)**c%m4
      f = c%A * 10.0_wp**lg
   end function zda_gofa

   ! Parse the ZDA INI-style config. Fills comps(1:n_comp).
   ! Optional ok (0-arg absent -> stop as before; present -> .false. on error).
   subroutine read_zda_config(path, n_comp, comps, ok)
      character(len=*),  intent(in)  :: path
      integer,           intent(out) :: n_comp
      type(zda_comp_t),  intent(out) :: comps(ZDA_MAXCOMP)
      logical, optional, intent(out) :: ok
      integer :: u, ios, ic, ieq, ihash
      character(len=256) :: line, key, val, prevcom
      character(len=256) :: raw

      if (present(ok)) ok = .true.
      n_comp = 0;  ic = 0;  prevcom = ''
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(ok)) then
            ok = .false.;  return
         else
            write(*,'(a,a)') ' read_zda_config: cannot open ', trim(path); stop 1
         end if
      end if
      do
         read(u, '(a)', iostat=ios) raw
         if (ios /= 0) exit
         line = adjustl(raw)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') then
            prevcom = trim(adjustl(line(2:)))     ! remember last comment (component label)
            cycle
         end if
         if (line(1:1) == '[') then
            if (index(line, '[Component') > 0) then
               ic = ic + 1
               if (ic > ZDA_MAXCOMP) then
                  if (present(ok)) then
                     close(u);  ok = .false.;  return
                  else
                     write(*,'(a)') ' read_zda_config: too many components'; stop 1
                  end if
               end if
               comps(ic)%label = trim(prevcom)    ! label comes from the line before? set on next #
            end if
            cycle
         end if
         ihash = index(line, '#');  if (ihash > 0) line = line(:ihash-1)
         ieq = index(line, '=');    if (ieq == 0) cycle
         key = trim(adjustl(line(:ieq-1)))
         val = trim(adjustl(line(ieq+1:)))
         if (trim(key) == 'Number of Components') then
            read(val, *) n_comp;  cycle
         end if
         if (ic == 0) cycle
         call set_key(comps(ic), key, val)
      end do
      close(u)
      if (n_comp == 0) n_comp = ic
   end subroutine read_zda_config

   subroutine set_key(c, key, val)
      type(zda_comp_t), intent(inout) :: c
      character(len=*), intent(in)    :: key, val
      select case (trim(key))
      case ('Cross Sections');  c%xsec  = trim(val)
      case ('Calorimetry');     c%calor = trim(val)
      case ('A');     read(val,*) c%A
      case ('a_min'); read(val,*) c%a_min
      case ('a_max'); read(val,*) c%a_max
      case ('c0');    read(val,*) c%c0
      case ('b0');    read(val,*) c%b0
      case ('b1');    read(val,*) c%b1
      case ('a1');    read(val,*) c%a1
      case ('m1');    read(val,*) c%m1
      case ('b2');    read(val,*) c%b2
      case ('a2');    read(val,*) c%a2
      case ('m2');    read(val,*) c%m2
      case ('b3');    read(val,*) c%b3
      case ('a3');    read(val,*) c%a3
      case ('m3');    read(val,*) c%m3
      case ('b4');    read(val,*) c%b4
      case ('a4');    read(val,*) c%a4
      case ('m4');    read(val,*) c%m4
      end select
   end subroutine set_key


   ! ------------------------------------------------------------------
   ! DustEM/Zubko Q-table reader. Header (NSIZE, NWAVE, density), then
   ! one block for each radius: "<a> = radius (micron)", a column-header line, and
   ! NWAVE rows of  x  lambda[um]  Q_abs  Q_sca  Q_ext  g.
   ! Returns a_um(nsize), lam_um(nwave) [um], qabs/qsca(nwave,nsize), rho[g/cm^3].
   ! ------------------------------------------------------------------
   subroutine read_zubko_optics(path, nsize, nwave, a_um, lam_um, qabs, qsca, rho, ok)
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: nsize, nwave
      real(wp), allocatable, intent(out) :: a_um(:), lam_um(:), qabs(:,:), qsca(:,:)
      real(wp),              intent(out) :: rho
      logical, optional,     intent(out) :: ok
      integer :: u, ios, ja, jw, idum
      real(wp) :: x, ldum
      logical  :: found
      character(len=256) :: line

      if (present(ok)) ok = .true.
      nsize = 0;  nwave = 0;  rho = 0.0_wp
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(ok)) then
            ok = .false.;  return
         else
            write(*,'(a,a)') ' read_zubko_optics: cannot open ', trim(path); stop 1
         end if
      end if

      ! --- header: scan until the first "= radius" line ---
      do
         read(u,'(a)', iostat=ios) line;  if (ios /= 0) exit
         if (index(line,'NSIZE') > 0) read(line(index(line,'NSIZE')+5:), *) nsize
         if (index(line,'NWAVE') > 0) read(line(index(line,'NWAVE')+5:), *) nwave
         if (index(line,'Density') > 0 .or. index(line,'gr/cm^3') > 0) &
            read(line, *) idum, rho
         if (index(line,'radius') > 0) exit          ! first radius block header
      end do
      if (nsize <= 0 .or. nwave <= 0) then
         if (present(ok)) then
            close(u);  ok = .false.;  return
         else
            write(*,'(a)') ' read_zubko_optics: failed to parse NSIZE/NWAVE'; stop 1
         end if
      end if
      allocate(a_um(nsize), lam_um(nwave), qabs(nwave,nsize), qsca(nwave,nsize))

      ! --- block 1 (line currently holds its "= radius" header) ---
      read(line, *) a_um(1)
      read(u,'(a)') line                              ! column header
      do jw = 1, nwave
         read(u, *) x, lam_um(jw), qabs(jw,1), qsca(jw,1)
      end do
      ! --- blocks 2..nsize ---
      do ja = 2, nsize
         call skip_to_radius(u, line, found)
         if (.not. found) then
            if (present(ok)) then
               close(u);  deallocate(a_um, lam_um, qabs, qsca);  ok = .false.;  return
            else
               write(*,'(a)') ' read_zubko_optics: unexpected EOF seeking radius'; stop 1
            end if
         end if
         read(line, *) a_um(ja)
         read(u,'(a)') line                           ! column header
         do jw = 1, nwave
            read(u, *) x, ldum, qabs(jw,ja), qsca(jw,ja)
         end do
      end do
      close(u)
   end subroutine read_zubko_optics

   ! Advance to the next "= radius" block header. found=.false. at EOF.
   subroutine skip_to_radius(u, line, found)
      integer,            intent(in)  :: u
      character(len=256), intent(out) :: line
      logical,            intent(out) :: found
      integer :: ios
      found = .false.
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) return
         if (index(line,'radius') > 0) then
            found = .true.;  return
         end if
      end do
   end subroutine skip_to_radius


   ! ------------------------------------------------------------------
   ! Specific enthalpy / heat-capacity calorimetry reader. Header:
   !   2 comment lines, then "MIN,TMAX,NT: <Tmin>, <Tmax> <NT>", then NT rows
   !   of  T[K]  U_spec[erg/gm]  C_spec[erg/gm/K].
   ! ------------------------------------------------------------------
   subroutine read_zubko_calor(path, nt, T, U_spec, C_spec, ok)
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: nt
      real(wp), allocatable, intent(out) :: T(:), U_spec(:), C_spec(:)
      logical, optional,     intent(out) :: ok
      integer :: u, ios, i, icol
      real(wp) :: tmin, tmax
      character(len=256) :: line

      if (present(ok)) ok = .true.
      nt = 0
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(ok)) then
            ok = .false.;  return
         else
            write(*,'(a,a)') ' read_zubko_calor: cannot open ', trim(path); stop 1
         end if
      end if
      do
         read(u,'(a)', iostat=ios) line;  if (ios /= 0) exit
         icol = index(line, ':')
         if (index(line,'NT') > 0 .and. icol > 0) then
            read(line(icol+1:), *) tmin, tmax, nt
            exit
         end if
      end do
      if (nt <= 0) then
         if (present(ok)) then
            close(u);  ok = .false.;  return
         else
            write(*,'(a)') ' read_zubko_calor: failed to parse NT'; stop 1
         end if
      end if
      allocate(T(nt), U_spec(nt), C_spec(nt))
      do i = 1, nt
         read(u, *) T(i), U_spec(i), C_spec(i)
      end do
      close(u)
   end subroutine read_zubko_calor


   ! ------------------------------------------------------------------
   ! Generic 2-column size-distribution table reader: a[um]  f(a).
   ! Skips any header line that does not start with a digit. f is whatever
   ! the file holds (e.g. dn/da or f(a)); the caller applies the bin/unit
   ! convention. Two passes: count rows, then read.
   ! ------------------------------------------------------------------
   subroutine read_dnda_table(path, n, a_um, fval, ok)
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: n
      real(wp), allocatable, intent(out) :: a_um(:), fval(:)
      logical, optional,     intent(out) :: ok
      integer :: u, ios, i
      real(wp) :: a, f
      character(len=256) :: line

      if (present(ok)) ok = .true.
      ! pass 1: count
      n = 0
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         if (present(ok)) then
            ok = .false.;  return
         else
            write(*,'(a,a)') ' read_dnda_table: cannot open ', trim(path); stop 1
         end if
      end if
      do
         read(u,'(a)', iostat=ios) line;  if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (.not. (line(1:1) >= '0' .and. line(1:1) <= '9')) cycle
         read(line, *, iostat=ios) a, f;  if (ios /= 0) cycle
         n = n + 1
      end do
      ! pass 2: read
      rewind(u)
      allocate(a_um(n), fval(n))
      i = 0
      do
         read(u,'(a)', iostat=ios) line;  if (ios /= 0) exit
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (.not. (line(1:1) >= '0' .and. line(1:1) <= '9')) cycle
         read(line, *, iostat=ios) a, f;  if (ios /= 0) cycle
         i = i + 1;  a_um(i) = a;  fval(i) = f
      end do
      close(u)
   end subroutine read_dnda_table

end module zubko_io
