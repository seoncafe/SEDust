module dustem_io
   ! Readers for the DustEM input-file formats, and the DustEM size
   ! distribution, for the two published models this tree carries as DustEM
   ! definitions:
   !
   !   data/themis/  THEMIS -- Jones et al. (2013), A&A 558, A62 (carbon)
   !                        + Koehler et al. (2014), A&A 565, L9 (silicates)
   !   data/g18d/    Guillet et al. (2018), A&A 610, A16, their Model D
   !
   ! Both are the models Hensley & Draine (2023), ApJ 948, 55, Sec. 6.2.2 and
   ! Fig. 16 compare their astrodust model against.  DustEM itself is
   ! Compiegne et al. (2011), A&A 525, A103.
   !
   ! Four file formats, one reader each:
   !
   !   GRAIN_*.DAT   the model definition: run keywords, the radiation-field
   !                 scaling G0, then one line per grain population giving its
   !                 size-distribution shape, mass per H, bulk density and
   !                 radius range.  read_dustem_grain.
   !   LAMBDA.DAT    the common wavelength grid every Q_/G_ table is written
   !                 against.  read_dustem_wavelengths.
   !   Q_*.DAT       Qabs and Qsca on (wavelength, radius).
   !                 read_dustem_qtable.
   !   G_*.DAT       the scattering asymmetry <cos theta> on the same grid.
   !                 read_dustem_gtable.
   !   C_*.DAT       the heat capacity per unit volume of grain material,
   !                 log10 C [erg K^-1 cm^-3] on (log10 T, radius).
   !                 read_dustem_heat_capacity.
   !
   ! Every one of them may carry '#' comment lines anywhere before a data
   ! block, so each reader steps over them and then reads free-format.
   !
   ! Two more entry points turn the tables into what a model needs on ITS own
   ! radius grid, both reproducing what DustEM does with the same tables:
   ! optics_at_radii (linear interpolation in radius) and
   ! grain_enthalpy_from_heat_capacity (the primitive of C dT).
   !
   ! The size-distribution shapes implemented here are the two DustEM writes
   ! analytically -- LOGN and PLAW, with the optional -ED exponential decay and
   ! -CV curvature factor on PLAW.  Anything else, the SIZE keyword included
   ! (which asks for a tabulated SIZE_*.DAT file), is refused by name rather
   ! than silently approximated.
   use constants, only: wp, pi
   use dust_model_mod, only: CHANNEL_NAME_LEN
   implicit none
   private
   public :: dustem_pop_t, DUSTEM_MAXPOP, DUSTEM_MAXPAR, DUSTEM_PROTON_MASS
   public :: read_dustem_grain, dustem_size_distribution, dustem_population_name
   public :: read_dustem_wavelengths, read_dustem_qtable, read_dustem_gtable
   public :: read_dustem_heat_capacity
   public :: optics_at_radii, grain_enthalpy_from_heat_capacity

   integer, parameter :: DUSTEM_MAXPOP = 16
   ! LOGN takes 2 shape parameters, PLAW 1, -ED 3 more and -CV 3 more; the
   ! shipped GRAIN files pad the line with trailing zeros beyond those, and the
   ! padding is read past, not used.
   integer, parameter :: DUSTEM_MAXPAR = 8

   ! Proton mass [g] as DustEM writes it (its `xmp` in DM_constants.f90).  It
   ! is the constant that turns the mass per H of a GRAIN_*.DAT line into a
   ! number per H, so reproducing a DustEM product means using DustEM's value
   ! and not a more recent CODATA one; the two differ by ~1e-8 relative.
   real(wp), parameter :: DUSTEM_PROTON_MASS = 1.67262158e-24_wp

   ! One grain population, i.e. one line of a GRAIN_*.DAT file.
   type :: dustem_pop_t
      character(len=64) :: gtype   = ''       ! names the Q_/G_/C_ files
      character(len=64) :: keyword = ''       ! the size-distribution keyword, lowercased
      integer  :: nsize = 0
      real(wp) :: mprop = 0.0_wp              ! M_dust/M_H of this population
      real(wp) :: rho   = 0.0_wp              ! bulk density [g/cm^3]
      real(wp) :: amin  = 0.0_wp, amax = 0.0_wp   ! [cm]
      real(wp) :: par(DUSTEM_MAXPAR) = 0.0_wp     ! shape parameters, in file order
      integer  :: npar  = 0                   ! how many of them this keyword uses
      logical  :: lognormal = .false.         ! LOGN
      logical  :: power_law = .false.         ! PLAW
      logical  :: exp_decay = .false.         ! -ED
      logical  :: curvature = .false.         ! -CV
      logical  :: polarized = .false.         ! -POL (no extra shape parameters)
   end type dustem_pop_t

contains

   ! ------------------------------------------------------------------
   ! What one population of a DustEM model is called when it has to be named
   ! on its own -- the /qtable group of the model's HDF5 product, written by
   ! sed/calc_qtable.x and read back by build_dustem.
   !
   ! The gtype is the natural name, and for every population of every shipped
   ! model but one it IS the name.  THEMIS is the exception: GRAIN_J13.DAT
   ! carries CM20 twice, once as the power law with exponential decay (the
   ! large a-C:H/a-C grains) and once as the log-normal (the small a-C
   ! grains), and the two have different radius ranges, so their stored optics
   ! are different arrays and cannot share a name.  A gtype that occurs more
   ! than once in the model therefore carries its size-distribution keyword,
   ! which is what distinguishes those two populations in the file itself.
   !
   ! The name is a function of the GRAIN_*.DAT alone, so the program that
   ! writes the product and the one that reads it agree without either
   ! recording a list.
   ! ------------------------------------------------------------------
   pure function dustem_population_name(pops, npop, ip) result(nm)
      type(dustem_pop_t), intent(in) :: pops(:)
      integer,            intent(in) :: npop, ip
      character(len=CHANNEL_NAME_LEN) :: nm
      integer :: j, nsame
      nsame = 0
      do j = 1, npop
         if (trim(pops(j)%gtype) == trim(pops(ip)%gtype)) nsame = nsame + 1
      end do
      if (nsame > 1) then
         nm = trim(pops(ip)%gtype)//'_'//trim(pops(ip)%keyword)
      else
         nm = trim(pops(ip)%gtype)
      end if
   end function dustem_population_name


   ! ------------------------------------------------------------------
   ! Position an open file at the next line that is neither blank nor a
   ! '#' comment, leaving it unread.  ios /= 0 at end of file.
   ! ------------------------------------------------------------------
   subroutine skip_comment_lines(u, ios)
      integer, intent(in)  :: u
      integer, intent(out) :: ios
      character(len=256)   :: line
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) return
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         backspace(u)
         return
      end do
   end subroutine skip_comment_lines


   ! ------------------------------------------------------------------
   ! Split a line on blanks and tabs.  nw = 0 for a blank line.
   ! ------------------------------------------------------------------
   subroutine split_words(line, maxw, nw, words)
      character(len=*), intent(in)  :: line
      integer,          intent(in)  :: maxw
      integer,          intent(out) :: nw
      character(len=*), intent(out) :: words(maxw)
      integer :: i, n, i0
      character(len=1), parameter :: TAB = achar(9)
      logical :: inword
      n = len_trim(line)
      nw = 0;  inword = .false.;  i0 = 1
      do i = 1, n
         if (line(i:i) == ' ' .or. line(i:i) == TAB) then
            if (inword) then
               if (nw < maxw) words(nw) = line(i0:i-1)
               inword = .false.
            end if
         else
            if (.not. inword) then
               nw = nw + 1;  i0 = i;  inword = .true.
            end if
         end if
      end do
      if (inword .and. nw <= maxw) words(nw) = line(i0:n)
   end subroutine split_words


   pure function lowercase(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s)) :: t
      integer :: i, k
      t = s
      do i = 1, len(s)
         k = iachar(s(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) t(i:i) = achar(k + 32)
      end do
   end function lowercase


   ! ------------------------------------------------------------------
   ! Trapezoid over a non-uniform abscissa -- DustEM's XINTEG2.  It is the
   ! quadrature the mass normalization below is defined by, so it is written
   ! out here rather than taken from elsewhere.
   ! ------------------------------------------------------------------
   pure function trapezoid_integral(x, y) result(s)
      real(wp), intent(in) :: x(:), y(:)
      real(wp) :: s
      integer  :: i
      s = 0.0_wp
      do i = 2, size(x)
         s = s + (x(i) - x(i-1)) * (y(i) + y(i-1)) * 0.5_wp
      end do
   end function trapezoid_integral


   ! ------------------------------------------------------------------
   ! GRAIN_*.DAT -- the model definition.
   !
   ! After the '#' comment lines: one line of run keywords (sdist, tempf, ...),
   ! which select DustEM's own outputs and mean nothing to a scalar model here;
   ! one line with G0, the radiation-field scaling; then one line per
   ! population,
   !
   !   <gtype> <nsize> <keyword> <Mdust/MH> <rho> <amin> <amax> <shape par...>
   !
   ! with amin, amax and the size parameters in CM and rho in g cm^-3.
   !
   ! ok present -> a bad file is reported through it; absent -> stop.
   ! ------------------------------------------------------------------
   subroutine read_dustem_grain(path, G0, npop, pops, ok)
      character(len=*),   intent(in)  :: path
      real(wp),           intent(out) :: G0
      integer,            intent(out) :: npop
      type(dustem_pop_t), intent(out) :: pops(DUSTEM_MAXPOP)
      logical, optional,  intent(out) :: ok
      integer, parameter :: MAXW = 32
      character(len=64)  :: words(MAXW)
      character(len=1024) :: line
      integer :: u, ios, nw, i, ip
      logical :: kok

      if (present(ok)) ok = .true.
      G0 = 0.0_wp;  npop = 0

      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call complain(' read_dustem_grain: cannot open '//trim(path), ok);  return
      end if

      ! line 1: run keywords -- read past.  They select which DustEM products a
      ! DustEM run writes; this tree builds the model, not those products.
      call skip_comment_lines(u, ios)
      if (ios /= 0) then
         close(u)
         call complain(' read_dustem_grain: no run-keyword line in '//trim(path), ok);  return
      end if
      read(u, '(a)', iostat=ios) line

      ! line 2: G0.  Read so that the file is fully accounted for; the scalar
      ! model is built from cross sections and a size distribution and never
      ! applies it.
      call skip_comment_lines(u, ios)
      if (ios /= 0) then
         close(u)
         call complain(' read_dustem_grain: no G0 line in '//trim(path), ok);  return
      end if
      read(u, *, iostat=ios) G0
      if (ios /= 0) then
         close(u)
         call complain(' read_dustem_grain: cannot read G0 from '//trim(path), ok);  return
      end if

      ! the population lines
      do
         call skip_comment_lines(u, ios)
         if (ios /= 0) exit
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         call split_words(line, MAXW, nw, words)
         if (nw < 7) then
            close(u)
            call complain(' read_dustem_grain: population line has fewer than 7 fields in ' &
                          //trim(path), ok)
            return
         end if
         if (npop >= DUSTEM_MAXPOP) then
            close(u)
            call complain(' read_dustem_grain: more populations than DUSTEM_MAXPOP in ' &
                          //trim(path), ok)
            return
         end if
         npop = npop + 1
         ip = npop
         pops(ip)%gtype   = trim(words(1))
         read(words(2), *, iostat=ios) pops(ip)%nsize
         if (ios /= 0) then
            close(u)
            call complain(' read_dustem_grain: bad nsize for '//trim(pops(ip)%gtype), ok)
            return
         end if
         pops(ip)%keyword = lowercase(trim(words(3)))
         read(words(4), *, iostat=ios) pops(ip)%mprop
         if (ios == 0) read(words(5), *, iostat=ios) pops(ip)%rho
         if (ios == 0) read(words(6), *, iostat=ios) pops(ip)%amin
         if (ios == 0) read(words(7), *, iostat=ios) pops(ip)%amax
         if (ios /= 0) then
            close(u)
            call complain(' read_dustem_grain: bad numeric field for ' &
                          //trim(pops(ip)%gtype), ok)
            return
         end if
         call parse_size_keyword(pops(ip), kok)
         if (.not. kok) then
            close(u)
            call complain(' read_dustem_grain: size-distribution keyword "' &
                          //trim(pops(ip)%keyword)//'" of grain type ' &
                          //trim(pops(ip)%gtype)//' is not implemented here.' &
                          //'  Only logn and plaw, with the optional -ed and -cv' &
                          //' factors on plaw, are.', ok)
            return
         end if
         if (nw < 7 + pops(ip)%npar) then
            close(u)
            call complain(' read_dustem_grain: too few shape parameters for ' &
                          //trim(pops(ip)%gtype), ok)
            return
         end if
         do i = 1, pops(ip)%npar
            read(words(7+i), *, iostat=ios) pops(ip)%par(i)
            if (ios /= 0) then
               close(u)
               call complain(' read_dustem_grain: bad shape parameter for ' &
                             //trim(pops(ip)%gtype), ok)
               return
            end if
         end do
         ! Any further field on the line is the zero padding the shipped G18D
         ! definition carries; DustEM reads only the parameters its keyword
         ! implies, and so does this.
      end do
      close(u)

      if (npop == 0) &
         call complain(' read_dustem_grain: no population lines in '//trim(path), ok)
   end subroutine read_dustem_grain


   ! Decode the dash-joined size-distribution keyword and set how many shape
   ! parameters it takes.  Returns kok = .false. for anything not implemented.
   subroutine parse_size_keyword(p, kok)
      type(dustem_pop_t), intent(inout) :: p
      logical,            intent(out)   :: kok
      character(len=64) :: rest, part
      integer :: k, nbase
      kok = .false.
      p%lognormal = .false.;  p%power_law = .false.
      p%exp_decay = .false.;  p%curvature = .false.;  p%polarized = .false.
      p%npar = 0

      rest = p%keyword
      ! base shape = the first dash-separated component
      k = index(rest, '-')
      if (k > 0) then
         part = rest(1:k-1);  rest = rest(k+1:)
      else
         part = rest;  rest = ''
      end if
      select case (trim(part))
      case ('logn');  p%lognormal = .true.;  nbase = 2
      case ('plaw');  p%power_law = .true.;  nbase = 1
      case default
         return                       ! 'size' and everything else: refused
      end select

      do
         if (len_trim(rest) == 0) exit
         k = index(rest, '-')
         if (k > 0) then
            part = rest(1:k-1);  rest = rest(k+1:)
         else
            part = rest;  rest = ''
         end if
         select case (trim(part))
         case ('ed');  p%exp_decay = .true.
         case ('cv');  p%curvature = .true.
         case ('pol'); p%polarized = .true.
         case default
            return                    ! chrg, zm, spin, dtls, beta, ...: refused
         end select
      end do

      ! DustEM applies -ED and -CV inside its PLAW branch only: its LOGN branch
      ! reads the keyword and then ignores both.  Refusing the combination here
      ! is deliberate -- accepting it would have to mean silently ignoring it.
      if (p%lognormal .and. (p%exp_decay .or. p%curvature)) return

      p%npar = nbase
      if (p%exp_decay) p%npar = p%npar + 3      ! at, ac, gamma
      if (p%curvature) p%npar = p%npar + 3      ! au, zeta, eta -- after the -ED three
      if (p%npar > DUSTEM_MAXPAR) return
      kok = .true.
   end subroutine parse_size_keyword


   ! ------------------------------------------------------------------
   ! The DustEM size distribution of one population.
   !
   ! The radius grid is nsize points spaced evenly in ln a between amin and
   ! amax, with the last point set to amax exactly so that rounding cannot put
   ! it outside the range the optics tables cover.
   !
   ! The shape function ava = a^4 dn/da is built in arbitrary units and then
   ! normalized so that the population carries the mass per H its GRAIN line
   ! states:
   !
   !   dn/da [cm^-1 per H] = ava / a^4 * m_p * (Mdust/MH) / masstot
   !   masstot             = trapezoid over ln a of  ava * rho * 4 pi / 3
   !
   ! With the trapezoid weights of the caller's size sum this makes
   ! sum_a (4/3) pi a^3 rho dn(a) come out at m_p * (Mdust/MH) exactly, which
   ! is what makes the model's dust mass per H the number its definition
   ! states.
   !
   ! a_cm, lna and dnda are allocated here to nsize.
   ! ------------------------------------------------------------------
   subroutine dustem_size_distribution(p, a_cm, lna, dnda, ok)
      type(dustem_pop_t),    intent(in)  :: p
      real(wp), allocatable, intent(out) :: a_cm(:), lna(:), dnda(:)
      real(wp), allocatable              :: ava(:)
      logical, optional,     intent(out) :: ok
      real(wp) :: dlna, arg, masstot, at, ac, gam, au, zeta, eta, zsign
      integer  :: n, j, k

      if (present(ok)) ok = .true.
      n = p%nsize
      if (n < 1 .or. p%amin <= 0.0_wp .or. p%amax <= 0.0_wp) then
         call complain(' dustem_size_distribution: bad size range for ' &
                       //trim(p%gtype), ok)
         return
      end if
      allocate(a_cm(n), lna(n), dnda(n), ava(n))

      if (n /= 1) then
         dlna = (log(p%amax) - log(p%amin)) / real(n-1, wp)
      else
         dlna = log(p%amax) - log(p%amin)
      end if
      do j = 1, n
         lna(j) = log(p%amin) + real(j-1, wp) * dlna
      end do
      ! DustEM overwrites the last node with log(amax) so that the grid ends on
      ! the stated maximum rather than on its rounded reconstruction.
      lna(n) = log(p%amax)
      a_cm   = exp(lna)
      a_cm(n) = p%amax

      ! --- shape function ava = a^4 dn/da, arbitrary units ---
      ! The -350 guard is DustEM's: it zeroes the term rather than letting the
      ! exponential underflow, so that the two codes agree bit for bit on which
      ! nodes are exactly zero.
      if (p%lognormal) then
         ! ava = exp( 3 ln a - (ln a - ln a0)^2 / (2 sigma^2) ).  The logarithm
         ! is NATURAL: DustEM's documentation says "log" but its code uses LOG.
         if (p%par(1) <= 0.0_wp .or. p%par(2) == 0.0_wp) then
            call complain(' dustem_size_distribution: lognormal centroid or sigma is 0 for ' &
                          //trim(p%gtype), ok)
            return
         end if
         do j = 1, n
            arg = 3.0_wp*lna(j) - 0.5_wp * ((lna(j) - log(p%par(1))) / p%par(2))**2
            if (arg > -350.0_wp) then
               ava(j) = exp(arg)
            else
               ava(j) = 0.0_wp
            end if
         end do
      else
         ! ava = a^(4 + alpha)
         do j = 1, n
            arg = (4.0_wp + p%par(1)) * lna(j)
            if (arg > -350.0_wp) then
               ava(j) = exp(arg)
            else
               ava(j) = 0.0_wp
            end if
         end do
         k = 1
         if (p%exp_decay) then
            ! exponential decay above a threshold radius: leave a < at alone
            at = p%par(k+1);  ac = p%par(k+2);  gam = p%par(k+3);  k = k + 3
            do j = 1, n
               if (a_cm(j) >= at) then
                  arg = -(((a_cm(j) - at) / ac)**gam)
                  if (arg > -350.0_wp) then
                     ava(j) = ava(j) * exp(arg)
                  else
                     ava(j) = 0.0_wp
                  end if
               end if
            end do
         end if
         if (p%curvature) then
            au = p%par(k+1);  zeta = abs(p%par(k+2));  eta = p%par(k+3)
            zsign = sign(1.0_wp, p%par(k+2))
            do j = 1, n
               ava(j) = ava(j) * (1.0_wp + zeta*(a_cm(j)/au)**eta)**zsign
            end do
         end if
      end if

      ! --- normalize to the stated mass per H ---
      masstot = trapezoid_integral(lna, ava * p%rho * 4.0_wp*pi/3.0_wp)
      if (masstot <= 0.0_wp) then
         call complain(' dustem_size_distribution: zero dust mass for '//trim(p%gtype), ok)
         return
      end if
      do j = 1, n
         dnda(j) = ava(j) / a_cm(j)**4 * DUSTEM_PROTON_MASS * p%mprop / masstot
      end do
      deallocate(ava)
   end subroutine dustem_size_distribution


   ! ------------------------------------------------------------------
   ! LAMBDA.DAT -- a count, then that many wavelengths in um, ascending,
   ! possibly several to a line.
   ! ------------------------------------------------------------------
   subroutine read_dustem_wavelengths(path, nwave, lam_um, ok)
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: nwave
      real(wp), allocatable, intent(out) :: lam_um(:)
      logical, optional,     intent(out) :: ok
      integer :: u, ios

      if (present(ok)) ok = .true.
      nwave = 0
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call complain(' read_dustem_wavelengths: cannot open '//trim(path), ok);  return
      end if
      call skip_comment_lines(u, ios)
      if (ios == 0) read(u, *, iostat=ios) nwave
      if (ios /= 0 .or. nwave < 2) then
         close(u)
         call complain(' read_dustem_wavelengths: cannot read the point count from ' &
                       //trim(path), ok)
         return
      end if
      allocate(lam_um(nwave))
      call skip_comment_lines(u, ios)
      if (ios == 0) read(u, *, iostat=ios) lam_um
      close(u)
      if (ios /= 0) then
         deallocate(lam_um);  nwave = 0
         call complain(' read_dustem_wavelengths: cannot read the grid from '//trim(path), ok)
      end if
   end subroutine read_dustem_wavelengths


   ! ------------------------------------------------------------------
   ! Q_*.DAT -- comments, nsize, nsize radii in um, then a QABS block of nwave
   ! rows of nsize values and a QSCA block of the same shape.  Row k of a block
   ! belongs to entry k of LAMBDA.DAT, so nwave comes from that file and is an
   ! input here.
   ! ------------------------------------------------------------------
   subroutine read_dustem_qtable(path, nwave, nsize, a_um, qabs, qsca, ok)
      character(len=*),      intent(in)  :: path
      integer,               intent(in)  :: nwave
      integer,               intent(out) :: nsize
      real(wp), allocatable, intent(out) :: a_um(:), qabs(:,:), qsca(:,:)
      logical, optional,     intent(out) :: ok
      integer :: u, ios

      if (present(ok)) ok = .true.
      nsize = 0
      call open_sized_table(path, u, nsize, a_um, ios)
      if (ios /= 0) then
         call complain(' read_dustem_qtable: cannot read the radii of '//trim(path), ok)
         return
      end if
      allocate(qabs(nwave, nsize), qsca(nwave, nsize))
      call read_lambda_block(u, nwave, nsize, qabs, ios)
      if (ios == 0) call read_lambda_block(u, nwave, nsize, qsca, ios)
      close(u)
      if (ios /= 0) then
         deallocate(a_um, qabs, qsca);  nsize = 0
         call complain(' read_dustem_qtable: cannot read the Q blocks of '//trim(path), ok)
      end if
   end subroutine read_dustem_qtable


   ! ------------------------------------------------------------------
   ! G_*.DAT -- the same layout with ONE block, the scattering asymmetry
   ! <cos theta>.
   ! ------------------------------------------------------------------
   subroutine read_dustem_gtable(path, nwave, nsize, a_um, gfac, ok)
      character(len=*),      intent(in)  :: path
      integer,               intent(in)  :: nwave
      integer,               intent(out) :: nsize
      real(wp), allocatable, intent(out) :: a_um(:), gfac(:,:)
      logical, optional,     intent(out) :: ok
      integer :: u, ios

      if (present(ok)) ok = .true.
      nsize = 0
      call open_sized_table(path, u, nsize, a_um, ios)
      if (ios /= 0) then
         call complain(' read_dustem_gtable: cannot read the radii of '//trim(path), ok)
         return
      end if
      allocate(gfac(nwave, nsize))
      call read_lambda_block(u, nwave, nsize, gfac, ios)
      close(u)
      if (ios /= 0) then
         deallocate(a_um, gfac);  nsize = 0
         call complain(' read_dustem_gtable: cannot read the g block of '//trim(path), ok)
      end if
   end subroutine read_dustem_gtable


   ! Open a Q_/G_ table and read its leading "nsize, then nsize radii" header.
   subroutine open_sized_table(path, u, nsize, a_um, ios)
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: u, nsize, ios
      real(wp), allocatable, intent(out) :: a_um(:)
      nsize = 0
      open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
      if (ios /= 0) return
      call skip_comment_lines(u, ios)
      if (ios == 0) read(u, *, iostat=ios) nsize
      if (ios /= 0 .or. nsize < 2) then
         close(u);  ios = 1;  return
      end if
      allocate(a_um(nsize))
      call skip_comment_lines(u, ios)
      if (ios == 0) read(u, *, iostat=ios) a_um
      if (ios /= 0) then
         deallocate(a_um);  close(u);  nsize = 0;  ios = 1
      end if
   end subroutine open_sized_table


   ! Read one block of nwave rows of nsize values, stepping over the '####'
   ! comment line that titles it.
   subroutine read_lambda_block(u, nwave, nsize, blk, ios)
      integer,  intent(in)  :: u, nwave, nsize
      real(wp), intent(out) :: blk(nwave, nsize)
      integer,  intent(out) :: ios
      integer :: k
      call skip_comment_lines(u, ios)
      if (ios /= 0) return
      do k = 1, nwave
         read(u, *, iostat=ios) blk(k, 1:nsize)
         if (ios /= 0) return
      end do
   end subroutine read_lambda_block


   ! ------------------------------------------------------------------
   ! C_*.DAT -- comments, nsize, nsize radii in um, ntemp, then ntemp rows of
   ! log10(T[K]) followed by log10(C[erg K^-1 cm^-3]) for each radius.  C is
   ! the heat capacity per unit volume of the grain material.
   ! ------------------------------------------------------------------
   subroutine read_dustem_heat_capacity(path, nsize, ntemp, a_um, log_T, log_C, ok)
      character(len=*),      intent(in)  :: path
      integer,               intent(out) :: nsize, ntemp
      real(wp), allocatable, intent(out) :: a_um(:), log_T(:), log_C(:,:)
      logical, optional,     intent(out) :: ok
      integer :: u, ios, i

      if (present(ok)) ok = .true.
      nsize = 0;  ntemp = 0
      call open_sized_table(path, u, nsize, a_um, ios)
      if (ios /= 0) then
         call complain(' read_dustem_heat_capacity: cannot read the radii of '//trim(path), ok)
         return
      end if
      call skip_comment_lines(u, ios)
      if (ios == 0) read(u, *, iostat=ios) ntemp
      ! The enthalpy normalization below anchors on tabulated points 2 and 6.
      if (ios /= 0 .or. ntemp < 6) then
         close(u);  deallocate(a_um);  nsize = 0
         call complain(' read_dustem_heat_capacity: needs at least 6 temperatures in ' &
                       //trim(path), ok)
         return
      end if
      allocate(log_T(ntemp), log_C(ntemp, nsize))
      call skip_comment_lines(u, ios)
      do i = 1, ntemp
         if (ios /= 0) exit
         read(u, *, iostat=ios) log_T(i), log_C(i, 1:nsize)
      end do
      close(u)
      if (ios /= 0) then
         deallocate(a_um, log_T, log_C);  nsize = 0;  ntemp = 0
         call complain(' read_dustem_heat_capacity: cannot read the C rows of '//trim(path), ok)
      end if
   end subroutine read_dustem_heat_capacity


   ! ------------------------------------------------------------------
   ! Put a tabulated optical efficiency -- Qabs, Qsca or <cos theta> -- on the
   ! model's own radius grid.
   !
   ! DustEM interpolates LINEARLY IN RADIUS at each wavelength independently
   ! (its INTPOL), so that is what happens here.  Outside the table DustEM
   ! returns zero; this refuses instead, because a model radius the optics do
   ! not cover is a broken model definition and a silent zero would hide it.
   ! Radii within rel_tol of an end of the table are pulled onto that end: the
   ! model grid is rebuilt from log(amin) and log(amax) and can land a few ulp
   ! outside a table written to the same values.
   ! ------------------------------------------------------------------
   subroutine optics_at_radii(a_tab_um, q_tab, a_um, q_out, ok, what)
      real(wp),         intent(in)  :: a_tab_um(:)      ! (ntab) ascending
      real(wp),         intent(in)  :: q_tab(:,:)       ! (nwave, ntab)
      real(wp),         intent(in)  :: a_um(:)          ! (na) model radii
      real(wp),         intent(out) :: q_out(:,:)       ! (nwave, na)
      logical, optional, intent(out) :: ok
      character(len=*), optional, intent(in) :: what    ! for the error message
      real(wp), parameter :: REL_TOL = 1.0e-8_wp
      real(wp) :: a, f
      integer  :: ntab, na, nwave, j, k, lo, hi, mid
      character(len=32) :: label

      if (present(ok)) ok = .true.
      label = 'optics';  if (present(what)) label = what
      ntab  = size(a_tab_um);  na = size(a_um);  nwave = size(q_tab, 1)

      do j = 1, na
         a = a_um(j)
         if (a < a_tab_um(1)) then
            if (a >= a_tab_um(1)*(1.0_wp - REL_TOL)) then
               a = a_tab_um(1)
            else
               call complain_radius(label, a, a_tab_um(1), a_tab_um(ntab), ok);  return
            end if
         else if (a > a_tab_um(ntab)) then
            if (a <= a_tab_um(ntab)*(1.0_wp + REL_TOL)) then
               a = a_tab_um(ntab)
            else
               call complain_radius(label, a, a_tab_um(1), a_tab_um(ntab), ok);  return
            end if
         end if
         lo = 1;  hi = ntab
         do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (a_tab_um(mid) <= a) then;  lo = mid;  else;  hi = mid;  end if
         end do
         f = (a - a_tab_um(lo)) / (a_tab_um(hi) - a_tab_um(lo))
         do k = 1, nwave
            q_out(k, j) = (1.0_wp - f)*q_tab(k, lo) + f*q_tab(k, hi)
         end do
      end do
   end subroutine optics_at_radii


   ! ------------------------------------------------------------------
   ! Total enthalpy U(T) of a grain of radius a, from the tabulated heat
   ! capacity per unit volume -- DustEM's GET_UC.
   !
   !   log10 C(i, a)  linear in RADIUS over the table radii, held constant
   !                  outside them (DustEM's INTPOL3)
   !   C(i)           = (4 pi / 3) a^3 * 10^log10C(i,a)          [erg/K]
   !   t(i)           = 10^logT(i)
   !   alpha          = dlogC/dlogT from tabulated points 2 and 6
   !   U(2)           = C(2) t(2) / (1 + alpha)
   !   U(i)           = U(2) + trapezoid from 2 to i of ln(10) C t d(log10 T)
   !
   ! i.e. the primitive of C dT written on the tabulated log10 T abscissa,
   ! anchored at the SECOND tabulated temperature, where the C ~ T^alpha limit
   ! of the low-temperature end fixes the constant of integration.
   !
   ! The result is then put on the caller's temperature grid by log-log
   ! interpolation, clamped to U(1) below the first tabulated temperature
   ! rather than extrapolated.
   ! ------------------------------------------------------------------
   subroutine grain_enthalpy_from_heat_capacity(a_um, a_tab_um, log_T, log_C, T_out, U_out, ok)
      real(wp),          intent(in)  :: a_um             ! grain radius [um]
      real(wp),          intent(in)  :: a_tab_um(:)      ! (ntab) ascending
      real(wp),          intent(in)  :: log_T(:)         ! (ntemp) log10 T [K]
      real(wp),          intent(in)  :: log_C(:,:)       ! (ntemp, ntab) log10 C [erg/K/cm^3]
      real(wp),          intent(in)  :: T_out(:)         ! (nt) target grid [K]
      real(wp),          intent(out) :: U_out(:)         ! (nt) enthalpy [erg]
      logical, optional, intent(out) :: ok
      real(wp), parameter :: REL_TOL = 1.0e-8_wp
      real(wp), allocatable :: lc(:), cc(:), tt(:), uu(:)
      real(wp) :: a, f, alpha, u0, acm3, lt, s, fr
      integer  :: ntab, ntemp, nt, i, j, lo, hi, mid

      if (present(ok)) ok = .true.
      ntab = size(a_tab_um);  ntemp = size(log_T);  nt = size(T_out)

      a = a_um
      if (a < a_tab_um(1)   .and. a >= a_tab_um(1)*(1.0_wp - REL_TOL))    a = a_tab_um(1)
      if (a > a_tab_um(ntab).and. a <= a_tab_um(ntab)*(1.0_wp + REL_TOL)) a = a_tab_um(ntab)
      if (a < a_tab_um(1) .or. a > a_tab_um(ntab)) then
         call complain_radius('heat capacity', a, a_tab_um(1), a_tab_um(ntab), ok);  return
      end if

      allocate(lc(ntemp), cc(ntemp), tt(ntemp), uu(ntemp))

      ! log10 C at this radius, linear in radius
      lo = 1;  hi = ntab
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (a_tab_um(mid) <= a) then;  lo = mid;  else;  hi = mid;  end if
      end do
      f = (a - a_tab_um(lo)) / (a_tab_um(hi) - a_tab_um(lo))
      do i = 1, ntemp
         lc(i) = (1.0_wp - f)*log_C(i, lo) + f*log_C(i, hi)
      end do

      acm3 = (4.0_wp*pi/3.0_wp) * (a * 1.0e-4_wp)**3        ! grain volume [cm^3]
      do i = 1, ntemp
         cc(i) = acm3 * 10.0_wp**lc(i)                      ! [erg/K]
         tt(i) = 10.0_wp**log_T(i)                          ! [K]
      end do

      alpha = (lc(6) - lc(2)) / (log_T(6) - log_T(2))
      u0    = cc(2) * tt(2) / (1.0_wp + alpha)

      ! cumulative trapezoid of ln(10) C T over log10 T, recentred on point 2
      uu = 0.0_wp
      do i = 2, ntemp
         uu(i) = uu(i-1) + (log_T(i) - log_T(i-1)) * 0.5_wp * log(10.0_wp) &
                 * (cc(i)*tt(i) + cc(i-1)*tt(i-1))
      end do
      s = uu(2)
      do i = 1, ntemp
         uu(i) = uu(i) - s + u0
      end do

      ! onto the caller's temperature grid, log-log
      do j = 1, nt
         if (T_out(j) <= tt(1)) then
            U_out(j) = uu(1)
            cycle
         end if
         lt = log(T_out(j))
         if (T_out(j) >= tt(ntemp)) then
            ! Above the table the last tabulated log-log slope is continued;
            ! the shipped C files reach 3000 K, past the sublimation of every
            ! material here, so this is a guard and not a working branch.
            fr = (log(uu(ntemp)) - log(uu(ntemp-1))) / (log(tt(ntemp)) - log(tt(ntemp-1)))
            U_out(j) = exp(log(uu(ntemp)) + fr*(lt - log(tt(ntemp))))
            cycle
         end if
         lo = 1;  hi = ntemp
         do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (tt(mid) <= T_out(j)) then;  lo = mid;  else;  hi = mid;  end if
         end do
         fr = (lt - log(tt(lo))) / (log(tt(hi)) - log(tt(lo)))
         U_out(j) = exp((1.0_wp - fr)*log(uu(lo)) + fr*log(uu(hi)))
      end do

      deallocate(lc, cc, tt, uu)
   end subroutine grain_enthalpy_from_heat_capacity


   ! Report an error through ok when the caller passed one, and stop otherwise
   ! -- the convention every reader in this tree follows.
   subroutine complain(msg, ok)
      character(len=*),  intent(in)  :: msg
      logical, optional, intent(out) :: ok
      if (present(ok)) then
         ok = .false.
         write(*,'(a)') trim(msg)
      else
         write(*,'(a)') trim(msg)
         stop 1
      end if
   end subroutine complain

   subroutine complain_radius(label, a, alo, ahi, ok)
      character(len=*),  intent(in)  :: label
      real(wp),          intent(in)  :: a, alo, ahi
      logical, optional, intent(out) :: ok
      character(len=200) :: s
      write(s,'(a,es13.6,a,es13.6,a,es13.6,a)') ' dustem_io: model radius ', a, &
         ' um is outside the '//trim(label)//' table, which covers ', alo, ' - ', ahi, ' um'
      call complain(trim(s), ok)
   end subroutine complain_radius

end module dustem_io
