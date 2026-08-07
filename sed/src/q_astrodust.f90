module q_astrodust_mod
   ! The DH21 astrodust refractive index, and Q(a, lambda) for the
   ! volume-equivalent sphere computed from it, for the extreme-ultraviolet
   ! band the T-matrix Q table does not cover.
   !
   ! MATERIAL.  index_DH21Ad_P0.20_0.00_1.400 is the Draine & Hensley (2021)
   ! astrodust effective dielectric function for porosity P = 0.20, iron
   ! fraction f_Fe = 0.00 and axial ratio b/a = 1.400.  It is the SAME file
   ! that tmatrix/driver/run_tmatrix.f90 reads to produce
   ! tmatrix/output/q_astrodust_P0.20_Fe0.00_1.400.dat, so this module and
   ! that table describe one material.  The Bruggeman mixing for (P, f_Fe)
   ! is already encoded in the file; nothing is mixed here.  The table spans
   ! E = 2.485e-5 to 1.2398e4 eV (lambda = 4.98883e4 to 1.00003e-4 um), so it
   ! covers the whole EUV.  set_astrodust_index_path selects a different
   ! (P, f_Fe, b/a); it must stay PAIRED with the Q table the caller reads,
   ! or the model silently changes material at the seam between the two.
   ! sed_init therefore takes the path as an argument and prints the file it
   ! is using whenever the EUV band is active.
   !
   ! WHAT THIS MODULE OFFERS.  Two separate things, and the EUV band chooses
   ! between them:
   !
   !   astrodust_index_at        the interpolated (n, k) at one wavelength and
   !                             nothing else.  This is the DEFAULT route: it
   !                             hands the index to the T-matrix
   !                             (spheroid_q of spheroid_optics), so the EUV band
   !                             is the same b/a = 1.400 oblate spheroid, on
   !                             the same random-orientation average, as the
   !                             Q table it continues.  Nothing about the
   !                             grain's shape changes at the seam.
   !
   !   q_astrodust_abs,          Bohren-Huffman Mie for the volume-equivalent
   !   q_astrodust_full          SPHERE of the same a_eff.  A deliberate shape
   !                             APPROXIMATION, kept because it costs
   !                             milliseconds where the T-matrix costs
   !                             minutes; sed_init selects it with
   !                             euv_tmatrix = .false.  Its domain of validity
   !                             is measured below.
   !
   ! SHAPE, AND THE DOMAIN OF VALIDITY OF THE MIE APPROXIMATION.  The Q table
   ! is a random-orientation T-matrix calculation for an OBLATE SPHEROID of
   ! axial ratio b/a = 1.400, indexed by the volume-equivalent radius a_eff.
   ! The Mie route feeds the same a_eff to the volume-equivalent SPHERE.  That
   ! substitution is valid only where it is used: at wavelengths SHORTER than
   ! the Q table's short-wavelength end (0.0912 um = 13.6 eV).
   !
   ! Measured against the T-matrix at the same (a, lambda), with
   ! dQabs = (Mie - T-matrix) / T-matrix:
   !
   !   a [um] \ E  13.6 eV   17.7     24.8     30       40       50       80      100
   !   4.73e-4      -0.37%   +0.71%   -0.01%   -0.08%   -0.13%   -0.02%   -0.01%  -0.00%
   !   1.0e-3       -0.20%   +0.75%   +0.05%   -0.00%   -0.03%   -0.02%   -0.01%  -0.01%
   !   5.0e-3       -0.23%   +0.75%   -0.09%   -0.17%   -0.24%   -0.23%   -0.10%  -0.08%
   !   1.0e-2       +0.26%   +0.24%   -0.60%   -0.57%   -0.37%   -0.27%   -0.22%  -0.20%
   !   5.0e-2       -1.28%   -1.47%   -1.59%   -1.58%   -1.47%   -1.27%   -0.97%  -0.77%
   !   1.0e-1       -1.62%   -1.74%   -1.85%   -1.87%   -1.86%   -1.77%   -1.55%   (GO)
   !   3.0e-1       -1.90%   -1.96%   -2.01%   -2.02%    (GO)     (GO)     (GO)    (GO)
   !
   ! ((GO) = the T-matrix run has handed over to its geometric-optics limit.)
   ! The error stays within ~2% over the whole band and the whole size range,
   ! but it does NOT fall monotonically with photon energy: only a <~ 2e-3 um
   ! improves monotonically, a >~ 0.01 um gets WORSE between 13.6 and 30 eV
   ! before recovering, and the largest sizes settle at a nearly
   ! energy-independent -1.9 to -2.0%.  The small and large limits are
   ! controlled by different physics, and neither of them is the "x >~ 1 with
   ! |m - 1| << 1" anomalous-diffraction argument:
   !
   !   * SMALL grains are in the RAYLEIGH limit, not the diffraction limit.
   !     The smallest size the astrodust distribution populates,
   !     a = 4.73e-4 um, has x = 2 pi a / lambda = 0.033 at the seam and still
   !     only 0.24 at 100 eV.  There the polarizability of an ellipsoid is
   !     alpha_j = V (eps - 1) / [4 pi (1 + L_j (eps - 1))], so the shape
   !     enters ONLY through the product L_j (eps - 1) of the depolarization
   !     factors with (eps - 1).  |eps - 1| falls from 1.84 at 13.6 eV to 1.10
   !     at 20 eV, 0.28 at 40 eV and 0.061 at 100 eV; as it does, the L_j drop
   !     out and sphere and spheroid converge on the same
   !     C_abs = (2 pi / lambda) V Im(eps).  That, not a large size parameter,
   !     is why these grains agree best and improve monotonically.
   !   * LARGE grains are in the geometric-optics limit, where C_ext -> 2 x
   !     (projected area) and the shape enters only through the
   !     orientation-averaged projected area.  By Cauchy's theorem that average
   !     is S/4, and for a b/a = 1.4 oblate spheroid S/4 exceeds the
   !     pi a_eff^2 of its volume-equivalent sphere by 2.12%.  Mie is thus low
   !     by 1 - 1/1.0212 = 2.08% in that limit -- which is exactly what the
   !     a >= 0.1 um rows converge to.  The error tends to a constant, NOT to
   !     zero, however far into the EUV one goes.
   !   * In between, the two trends cross; that is where the non-monotonic
   !     rows come from.
   !
   ! So: bounded by ~2%, smallest for the small grains that carry most of the
   ! EUV absorption, and asymptotically a fixed ~2% underestimate for the
   ! largest ones -- an error that does not vanish however far into the EUV
   ! the band is carried, which is why the T-matrix and not this is the
   ! default.  The approximation must NOT be used at longer wavelengths, where
   ! the resonances and the b/a = 1.4 depolarization factors make the spheroid
   ! and the sphere genuinely different; there the T-matrix table is the
   ! model's optics.
   !
   ! The eV -> um conversion uses hc = 1.23984 eV um, matching q_silicate_mod
   ! and q_graphite_mod.  run_tmatrix.f90's reader rounds the same constant to
   ! 1.2398, i.e. the two wavelength axes differ by 3e-5 in relative terms.
   ! On this module's own EUV grid over 13.6-100 eV that shift moves Q_abs by
   ! 2e-5 (median) and at most 3e-4; a wavelength landing inside one of the
   ! table's absorption-edge steps can reach 5e-3, because the tabulation
   ! itself steps there (the largest such step below 100 eV is 5% in k across
   ! 0.03% in E, at 66.0 eV).  The table's ENERGY spacing is 3.6% at 13.6 eV
   ! and 13% near 40 eV, with 0.02-0.6% patches bracketing the edges, so the
   ! edges and that spacing -- not the constant -- set the resolution here.
   !
   ! Structure (parsing, caching, interpolation, Mie call) mirrors
   ! q_silicate_mod, which does the same job for the D03 astrosilicate.  The
   ! one difference is that the table length is counted at load time rather
   ! than fixed as a parameter, because the file path is settable.
   !
   ! THREAD SAFETY.  load_astrodust_index WRITES the cached table; everything
   ! else only reads it.  A caller that evaluates the optics from several
   ! threads must therefore load the table before the threads start (sed_init
   ! does), after which any number of threads may call astrodust_index_at,
   ! q_astrodust_abs or q_astrodust_full at once.

   use constants,   only: wp
   use sed_mathlib, only: interp
   use mie_mod,     only: mie
   use sed_paths,   only: sed_data_path
   implicit none
   private
   public :: astrodust_index_at
   public :: q_astrodust_abs
   public :: q_astrodust_full
   public :: set_astrodust_index_path
   public :: get_astrodust_index_path
   public :: load_astrodust_index
   public :: astrodust_index_lambda_range

   character(len=*), parameter :: F_AD_DEFAULT = &
        'dielectric/index_DH21Ad_P0.20_0.00_1.400'
   real(wp), parameter :: PI_LOC = 3.141592653589793238462643383279502884197_wp

   ! Blank means "the default, wherever the data root currently points".  A
   ! path a caller named is stored here as given and never composed with the
   ! root, so an absolute path from a host stays absolute.
   character(len=512) :: f_ad   = ''
   logical  :: loaded = .false.
   integer  :: nad    = 0
   real(wp), allocatable :: ad_eV(:), ad_n(:), ad_k(:), ad_wavl(:)

contains

   subroutine set_astrodust_index_path(path)
      ! Select a different astrodust dielectric function (another porosity,
      ! iron fraction or axial ratio).  The cached table is dropped, so the
      ! next optics call reloads from the new path.
      character(len=*), intent(in) :: path

      if (trim(path) == trim(get_astrodust_index_path())) return
      f_ad = path
      loaded = .false.
      if (allocated(ad_eV)) deallocate(ad_eV, ad_n, ad_k, ad_wavl)
      nad = 0
   end subroutine set_astrodust_index_path


   function get_astrodust_index_path() result(path)
      ! The dielectric function currently selected, so a caller can report
      ! which material its EUV band is made of.  Resolved, not stored: the
      ! default follows the data root.
      character(len=512) :: path
      if (len_trim(f_ad) == 0) then
         path = sed_data_path(F_AD_DEFAULT)
      else
         path = f_ad
      end if
   end function get_astrodust_index_path


   subroutine astrodust_index_lambda_range(lam_lo, lam_hi)
      ! Shortest and longest wavelength [um] the loaded table covers.  Outside
      ! this range `interp` returns the boundary value, i.e. a CONSTANT (n, k),
      ! so a caller that needs optics beyond it must refuse rather than let the
      ! frozen index pass as physics.
      real(wp), intent(out) :: lam_lo, lam_hi

      if (.not. loaded) call load_astrodust_index()
      lam_lo = minval(ad_wavl)
      lam_hi = maxval(ad_wavl)
   end subroutine astrodust_index_lambda_range


   subroutine load_astrodust_index(ok)
      ! index_DH21Ad_* columns: E[eV]  Re(m)-1  Im(m)  Re(eps)-1  Im(eps)
      ! (2 header lines, rows in ascending E = descending lambda).
      ! n, k are read directly: n = 1 + col2, k = col3.
      !
      ! Optional ok: absent -> message + stop on any error (what the direct
      ! callers below rely on); present -> return .false. with the module left
      ! unloaded instead of stopping, so sed_init can report the failure
      ! through its `status` argument rather than killing an RT host.
      logical, optional, intent(out) :: ok
      integer  :: i, u, ios
      real(wp) :: ener, rn1, rk, e1, e2
      character(len=512) :: line

      if (present(ok)) ok = .true.
      if (loaded) return

      open(newunit=u, file=get_astrodust_index_path(), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a,a)') 'q_astrodust: cannot open ', trim(get_astrodust_index_path())
         if (present(ok)) then
            ok = .false.;  return
         else
            stop 1
         end if
      end if

      ! Count the data rows, then rewind and read them.
      read(u, '(a)') line
      read(u, '(a)') line
      nad = 0
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         nad = nad + 1
      end do
      if (nad < 2) then
         write(*,'(a,i0,a,a)') 'q_astrodust: only ', nad, ' data rows in ', trim(get_astrodust_index_path())
         close(u)
         nad = 0
         if (present(ok)) then
            ok = .false.;  return
         else
            stop 1
         end if
      end if
      rewind(u)

      allocate(ad_eV(nad), ad_n(nad), ad_k(nad), ad_wavl(nad))
      read(u, '(/)')
      do i = 1, nad
         read(u, *, iostat=ios) ener, rn1, rk, e1, e2
         if (ios /= 0) then
            write(*,'(a,i0,a,a)') 'q_astrodust: malformed data row ', i, ' in ', trim(get_astrodust_index_path())
            close(u)
            deallocate(ad_eV, ad_n, ad_k, ad_wavl)
            nad = 0
            if (present(ok)) then
               ok = .false.;  return
            else
               stop 1
            end if
         end if
         ad_eV(i)   = ener
         ad_n(i)    = 1.0_wp + rn1
         ad_k(i)    = rk
         ad_wavl(i) = 1.23984_wp / ener     ! [um]
      end do
      close(u)
      loaded = .true.
   end subroutine load_astrodust_index


   subroutine astrodust_index_at(lambda, nr, ki)
      ! Refractive index m = nr + i*ki of the loaded astrodust dielectric
      ! function at one wavelength [um].  Exactly the interpolation
      ! q_astrodust_abs and q_astrodust_full do internally, with the Mie call
      ! left out, so a caller that brings its own optics -- the T-matrix EUV
      ! band of sed_astrodust -- reads the SAME index the sphere route would.
      !
      ! The lazy load below WRITES module state and so must never be reached
      ! from inside a parallel region: call load_astrodust_index() once before
      ! the threads start (sed_init does).  Past that point this routine only
      ! reads the cached table and is safe to call from many threads at once.
      real(wp), intent(in)  :: lambda
      real(wp), intent(out) :: nr, ki

      if (.not. loaded) call load_astrodust_index()

      ! ad_wavl is descending (eV ascending in file); interp handles both.
      call interp(ad_wavl, ad_n, lambda, nr)
      call interp(ad_wavl, ad_k, lambda, ki)
   end subroutine astrodust_index_at


   subroutine q_astrodust_abs(agrain, lambda, Qabs)
      ! Q_abs for an astrodust sphere of volume-equivalent radius agrain.
      ! agrain, lambda: um. Qabs: C_abs/(pi a^2).
      real(wp), intent(in)  :: agrain, lambda
      real(wp), intent(out) :: Qabs
      real(wp) :: x, nr, ki, Qext1, Qsca1, alb1, gsca1

      if (.not. loaded) call load_astrodust_index()

      ! ad_wavl is descending (eV ascending in file); interp handles both.
      call interp(ad_wavl, ad_n, lambda, nr)
      call interp(ad_wavl, ad_k, lambda, ki)

      x = 2.0_wp * PI_LOC * agrain / lambda
      call mie(nr, ki, x, Qext1, Qsca1, Qabs, alb1, gsca1)
   end subroutine q_astrodust_abs


   subroutine q_astrodust_full(agrain, lambda, Qext, Qsca, Qabs, gsca)
      ! Full Mie output (extinction, scattering, absorption efficiencies and
      ! scattering asymmetry g) for an astrodust sphere of volume-equivalent
      ! radius agrain. agrain, lambda: um. Same dielectric path as
      ! q_astrodust_abs -- just keeps every Mie return.  This is the sphere
      ! APPROXIMATION of the module header, reached only when a caller asks
      ! for it (sed_init's euv_tmatrix = .false.).
      real(wp), intent(in)  :: agrain, lambda
      real(wp), intent(out) :: Qext, Qsca, Qabs, gsca
      real(wp) :: x, nr, ki, alb1

      if (.not. loaded) call load_astrodust_index()

      call interp(ad_wavl, ad_n, lambda, nr)
      call interp(ad_wavl, ad_k, lambda, ki)

      x = 2.0_wp * PI_LOC * agrain / lambda
      call mie(nr, ki, x, Qext, Qsca, Qabs, alb1, gsca)
   end subroutine q_astrodust_full

end module q_astrodust_mod
