program calc_polext
   !====================================================================
   ! Polarized extinction cross section per H for the HD23 astrodust
   ! Milky-Way model, checked against the published release file.
   !
   !   p^max/N_H (lambda) = sum_a (dn_Ad/N_H)(a) * C_pol,ext(lambda,a)
   !                              * f_align(a)                  [cm^2/H]
   !
   ! with, from the orientation-resolved DH21 spheroid table,
   !   C_pol,ext = 0.5 * (Q_ext(jori=3) - Q_ext(jori=2)) * pi * a_eff^2
   ! i.e. the difference between E perp a and E || a at k perp a, halved,
   ! which is the extinction difference a perfectly aligned grain with its
   ! symmetry axis in the plane of the sky presents to the two linear
   ! polarizations.
   !
   ! Only astrodust grains are summed. HD23 take the PAHs to be
   ! unaligned (f_align = 0), so they contribute nothing here.
   !
   ! The size distribution dn_Ad/N_H comes from size_distribution.dat and
   ! is ALREADY integrated over each size bin -- do not multiply by da.
   !
   ! Output: output/polarized_extinction_ours.dat
   !   lambda[um]  ours[cm^2/H]  reference[cm^2/H]  ratio
   !
   ! The reference file goes negative below ~0.11 um, so the reported
   ! median/max deviations are taken over lambda >= 0.11 um only.
   !====================================================================
   use constants,        only: wp, pi, um2cm
   use q_table_jori_mod, only: load_q_table_jori, falign_hd23, &
                               nj_lam, nj_aeff, lam_j, aeff_j, qpol_ext
   use size_dist_mod,    only: load_size_dist, n_size, a_dist, dn_ad
   implicit none

   character(len=*), parameter :: F_Q    = '../data/dielectric/q_DH21Ad_P0.20_Fe0.00_1.400.dat.gz'
   character(len=*), parameter :: F_WAVE = '../data/dielectric/DH21_wave'
   character(len=*), parameter :: F_AEFF = '../data/dielectric/DH21_aeff'
   character(len=*), parameter :: F_SD   = '../data/release/size_distribution.dat'
   character(len=*), parameter :: F_REF  = '../data/release/polarized_extinction.dat'
   character(len=*), parameter :: F_OUT  = 'output/polarized_extinction_ours.dat'

   real(wp), parameter :: LAM_MIN_CMP = 0.11_wp   ! reference is negative below this

   real(wp), allocatable :: p_grid(:)        ! (nj_lam) on the DH21 wavelength grid
   real(wp), allocatable :: lam_ref(:), p_ref(:), p_ours(:)
   real(wp), allocatable :: dev(:)
   real(wp), allocatable :: wt(:)            ! (n_size) dn * f_align * pi a^2
   real(wp) :: qi, x, t, med, dmax, ratio
   integer  :: n_ref, jw, ia, i, lo, hi, mid, u, ios, ncmp
   logical  :: ok
   character(len=512) :: line

   write(*,'(a)') ' calc_polext: loading orientation-resolved DH21 optics ...'
   call load_q_table_jori(F_Q, F_WAVE, F_AEFF, ok)
   if (.not. ok) then
      write(*,'(a)') ' calc_polext: failed to load the DH21 Q table.'
      stop 1
   end if
   write(*,'(a,i0,a,i0)') '   NLAM=', nj_lam, '  NA=', nj_aeff

   call load_size_dist(F_SD, ok)
   if (.not. ok) then
      write(*,'(a)') ' calc_polext: failed to load the size distribution.'
      stop 1
   end if
   write(*,'(a,i0)') '   size bins=', n_size

   ! ---- size weights: dn/H * f_align * geometric cross section --------
   allocate(wt(n_size))
   do ia = 1, n_size
      wt(ia) = dn_ad(ia) * falign_hd23(a_dist(ia)) &
               * pi * (a_dist(ia) * um2cm)**2
   end do

   ! ---- size integral at each tabulated wavelength ---------------------
   ! Q_pol is interpolated log-linearly in a onto the size-distribution
   ! grid. That grid is the DH21 grid truncated at the small end and
   ! rounded to 5 digits, so the interpolation is very nearly the identity.
   allocate(p_grid(nj_lam))
   p_grid = 0.0_wp
   do jw = 1, nj_lam
      do ia = 1, n_size
         call interp_a(jw, a_dist(ia), qi)
         p_grid(jw) = p_grid(jw) + wt(ia) * qi
      end do
   end do

   ! ---- reference file --------------------------------------------------
   call count_data(F_REF, n_ref)
   if (n_ref < 2) then
      write(*,'(a)') ' calc_polext: cannot read '//F_REF
      stop 1
   end if
   allocate(lam_ref(n_ref), p_ref(n_ref), p_ours(n_ref))
   open(newunit=u, file=F_REF, status='old', action='read')
   i = 0
   do
      read(u,'(a)',iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#')    cycle
      i = i + 1
      read(line,*) lam_ref(i), p_ref(i)
   end do
   close(u)

   ! ---- put our result on the reference wavelengths (linear in log lam) --
   do i = 1, n_ref
      x = log(lam_ref(i))
      if (x <= log(lam_j(1))) then
         p_ours(i) = p_grid(1)
      else if (x >= log(lam_j(nj_lam))) then
         p_ours(i) = p_grid(nj_lam)
      else
         lo = 1;  hi = nj_lam
         do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (log(lam_j(mid)) <= x) then
               lo = mid
            else
               hi = mid
            end if
         end do
         t = (x - log(lam_j(lo))) / (log(lam_j(hi)) - log(lam_j(lo)))
         p_ours(i) = (1.0_wp - t) * p_grid(lo) + t * p_grid(hi)
      end if
   end do

   ! ---- write -----------------------------------------------------------
   open(newunit=u, file=F_OUT, status='replace', action='write')
   write(u,'(a)') '# Polarized extinction per H: this code vs HD23 release'
   write(u,'(a)') '# lambda[um]  ours[cm^2/H]  reference[cm^2/H]  ratio'
   do i = 1, n_ref
      if (p_ref(i) /= 0.0_wp) then
         ratio = p_ours(i) / p_ref(i)
      else
         ratio = 0.0_wp
      end if
      write(u,'(4es14.6)') lam_ref(i), p_ours(i), p_ref(i), ratio
   end do
   close(u)
   write(*,'(a)') ' calc_polext: wrote '//F_OUT

   ! ---- deviation statistics over the valid range ------------------------
   allocate(dev(n_ref))
   ncmp = 0
   do i = 1, n_ref
      if (lam_ref(i) < LAM_MIN_CMP) cycle
      if (p_ref(i) <= 0.0_wp)       cycle
      ncmp = ncmp + 1
      dev(ncmp) = abs(p_ours(i) / p_ref(i) - 1.0_wp)
   end do

   if (ncmp > 0) then
      call sort_head(dev, ncmp)
      if (mod(ncmp, 2) == 0) then
         med = 0.5_wp * (dev(ncmp/2) + dev(ncmp/2 + 1))
      else
         med = dev(ncmp/2 + 1)
      end if
      dmax = dev(ncmp)
      write(*,'(a,i0,a,f6.3,a)') ' calc_polext: ', ncmp, &
         ' points with lambda >= ', LAM_MIN_CMP, ' um'
      write(*,'(a,f8.4,a)') '   median |ours/ref - 1| = ', 100.0_wp*med,  ' %'
      write(*,'(a,f8.4,a)') '   max    |ours/ref - 1| = ', 100.0_wp*dmax, ' %'
   end if

   deallocate(wt, p_grid, lam_ref, p_ref, p_ours, dev)

contains

   subroutine interp_a(jw, a_target, q)
      ! Log-linear interpolation of qpol_ext in a_eff at fixed lambda index.
      integer,  intent(in)  :: jw
      real(wp), intent(in)  :: a_target
      real(wp), intent(out) :: q
      integer  :: l, h, m
      real(wp) :: xa, ta

      if (a_target <= aeff_j(1)) then
         q = qpol_ext(jw, 1);        return
      end if
      if (a_target >= aeff_j(nj_aeff)) then
         q = qpol_ext(jw, nj_aeff);  return
      end if
      xa = log(a_target)
      l = 1;  h = nj_aeff
      do while (h - l > 1)
         m = (l + h) / 2
         if (log(aeff_j(m)) <= xa) then
            l = m
         else
            h = m
         end if
      end do
      ta = (xa - log(aeff_j(l))) / (log(aeff_j(h)) - log(aeff_j(l)))
      q  = (1.0_wp - ta) * qpol_ext(jw, l) + ta * qpol_ext(jw, h)
   end subroutine interp_a


   subroutine count_data(filename, n)
      ! Count non-blank, non-`#` lines.
      character(len=*), intent(in)  :: filename
      integer,          intent(out) :: n
      integer :: uu, is
      character(len=512) :: ln
      n = 0
      open(newunit=uu, file=filename, status='old', action='read', iostat=is)
      if (is /= 0) return
      do
         read(uu,'(a)',iostat=is) ln
         if (is /= 0) exit
         ln = adjustl(ln)
         if (len_trim(ln) == 0) cycle
         if (ln(1:1) == '#')    cycle
         n = n + 1
      end do
      close(uu)
   end subroutine count_data


   subroutine sort_head(v, n)
      ! Ascending insertion sort of v(1:n); n is ~1e3 here.
      real(wp), intent(inout) :: v(:)
      integer,  intent(in)    :: n
      integer  :: ii, jj
      real(wp) :: key
      do ii = 2, n
         key = v(ii)
         jj  = ii - 1
         do while (jj >= 1)
            if (v(jj) <= key) exit
            v(jj+1) = v(jj)
            jj = jj - 1
         end do
         v(jj+1) = key
      end do
   end subroutine sort_head

end program calc_polext
