program test_dustem_polarized_extinction
   ! Compare the polarized extinction build_dustem forms for a DustEM-defined
   ! model against DustEM's own, for the very same input files:
   !
   !   G18 Model D   ../data/g18d/reference/EXT_POL_G17_ModelD.RES   see where.txt
   !
   ! The reference is a DustEM run on the GRAIN file, the Q1_/Q2_ tables and
   ! the ALIGN file this tree builds the model from, so the comparison isolates
   ! whether this tree reproduces DustEM's polarization conventions -- the
   ! 0.5*(Q2 - Q1) difference, the alignment law read off ALIGN_*.DAT, and the
   ! fact that the efficiency weights the size integral rather than the cross
   ! sections.  It never touches the solver: extinction is the size integral
   ! alone.
   !
   ! Reference format, the same as EXT_*.RES: '#' header, then a line
   ! "ntype nwave", then nwave rows of
   !   lambda[um]  ABS(1..ntype)  SCA(1..ntype)  PEXT(total)
   ! with every cross section in cm^2 per H multiplied by 1e21.  ABS and SCA
   ! are DustEM's EXTINCTION_POL, i.e. the size integrals of
   ! 0.5*(Q2_abs - Q1_abs)*f_pol and 0.5*(Q2_sca - Q1_sca)*f_pol, and PEXT is
   ! their sum over the populations.  Here the same three quantities are
   !
   !   ABS(ip)  = sum_a dn(a) f_align(a) Cpol(lam,a)
   !   SCA(ip)  = sum_a dn(a) f_align(a) [Cpol_ext(lam,a) - Cpol(lam,a)]
   !   PEXT     = size_integrated_extinction(..., Cpol_ext = ...)
   !
   ! The wavelengths are LAMBDA.DAT, so the model is built on the whole DustEM
   ! grid (include_euv) and no interpolation enters the comparison.
   !
   ! Two more things are checked here, because both are properties of the same
   ! build:
   !
   !   * THEMIS, which declares no aligned population, carries no polarized
   !     optics at all and its polarized extinction is identically zero.
   !   * dust_set_alignment overrides the tabulated DustEM efficiency with the
   !     HD23 power-law rolloff on every aligned population, as it does for
   !     astrodust, and clears align_tabulated.
   !
   ! Passes when every population's ABS and SCA and the total PEXT agree to
   ! 1e-4 relative -- far inside the reference file's own 7-digit print
   ! precision, which they reach.
   use constants,         only: wp
   use sed_astrodust_mod, only: dust_model_t, build_dustem, &
                                size_integrated_extinction, &
                                dust_has_polarized_optics
   use dust_model_mod,    only: dust_set_alignment
   use size_dist_mod,     only: falign_powerlaw
   implicit none

   integer,  parameter :: NT_IN  = 100
   real(wp), parameter :: T_LO = 1.0_wp, T_HI = 3000.0_wp
   real(wp), parameter :: TOL = 1.0e-4_wp
   real(wp), parameter :: SCALE = 1.0e21_wp

   character(len=*), parameter :: G18_GRAIN = '../data/g18d/GRAIN_G17_ModelD.DAT'
   character(len=*), parameter :: G18_DDIR  = '../data/g18d/'
   character(len=*), parameter :: G18_REF   = &
        '../data/g18d/reference/EXT_POL_G17_ModelD.RES'
   character(len=*), parameter :: THEMIS_GRAIN = '../data/themis/GRAIN_J13.DAT'
   character(len=*), parameter :: THEMIS_DDIR  = '../data/themis/'

   logical :: ok, all_ok

   all_ok = .true.

   call compare_against_dustem('G18 Model D', G18_GRAIN, G18_DDIR, G18_REF, ok)
   all_ok = all_ok .and. ok
   write(*,'(a)') ''

   call check_unaligned_model('THEMIS', THEMIS_GRAIN, THEMIS_DDIR, ok)
   all_ok = all_ok .and. ok
   write(*,'(a)') ''

   if (.not. all_ok) then
      write(*,'(a,es10.2)') ' FAIL: tolerance is ', TOL;  stop 1
   end if
   write(*,'(a)') ' PASS'

contains

   subroutine compare_against_dustem(label, grain_path, data_dir, ref_path, ok)
      character(len=*), intent(in)  :: label, grain_path, data_dir, ref_path
      logical,          intent(out) :: ok

      type(dust_model_t)    :: m
      real(wp), allocatable :: t(:,:), Cext(:), Cabs(:), Csca(:), Cpolext(:)
      real(wp), allocatable :: pabs_pop(:,:), psca_pop(:,:)
      real(wp) :: dev, w
      integer  :: ntype, nwave, ip, ja, jw, status, npop

      ok = .true.
      write(*,'(a)') ' === '//label//', polarized extinction'

      call read_dustem_res(ref_path, ntype, nwave, t, ok)
      if (.not. ok) return

      ! ---- the model, on the whole DustEM grid ---------------------------
      call build_dustem(m, grain_path, data_dir, NT_IN, T_LO, T_HI, status, &
                        include_euv=.true.)
      if (status /= 0) then
         write(*,'(a,i0)') '   build_dustem failed, status = ', status
         deallocate(t);  ok = .false.;  return
      end if
      npop = size(m%pops)
      write(*,'(a,i0,a,i0,a)') '   ', npop, ' populations, ', m%NLAM, ' wavelengths'

      if (.not. dust_has_polarized_optics(m)) then
         write(*,'(a)') '   the model carries no polarized optics'
         deallocate(t);  ok = .false.;  return
      end if
      if (.not. m%align_tabulated) then
         write(*,'(a)') '   align_tabulated is .false., but the alignment law' // &
              ' is not the HD23 power law'
         ok = .false.
      end if
      if (npop /= ntype) then
         write(*,'(a,i0,a,i0)') '   population count: model ', npop, ', reference ', ntype
         ok = .false.
      end if
      if (m%NLAM /= nwave) then
         write(*,'(a,i0,a,i0)') '   wavelength count: model ', m%NLAM, ', reference ', nwave
         ok = .false.
      end if
      if (.not. ok) then
         deallocate(t);  return
      end if

      ! ---- the size integrals --------------------------------------------
      allocate(Cext(nwave), Cabs(nwave), Csca(nwave), Cpolext(nwave))
      call size_integrated_extinction(m, Cext, Cabs, Csca, Cpol_ext=Cpolext)

      ! The absorption and scattering halves, population by population. The
      ! model stores the polarized ABSORPTION (Cpol) and the polarized
      ! EXTINCTION (Cpol_ext); their difference is the polarized scattering,
      ! which is the second block of the reference.
      allocate(pabs_pop(nwave, npop), psca_pop(nwave, npop))
      pabs_pop = 0.0_wp;  psca_pop = 0.0_wp
      do ip = 1, npop
         if (.not. allocated(m%pops(ip)%Cpol)) cycle
         do ja = 1, size(m%pops(ip)%dn)
            w = m%pops(ip)%dn(ja) * m%pops(ip)%falign(ja)
            do jw = 1, nwave
               pabs_pop(jw, ip) = pabs_pop(jw, ip) + w * m%pops(ip)%Cpol(jw, ja)
               psca_pop(jw, ip) = psca_pop(jw, ip) &
                    + w * (m%pops(ip)%Cpol_ext(jw, ja) - m%pops(ip)%Cpol(jw, ja))
            end do
         end do
      end do

      ! ---- compare --------------------------------------------------------
      write(*,'(a)') '   max |model/reference - 1|'
      call report_column('PEXT', '', Cpolext*SCALE, t(2*ntype+2,:), ok)
      do ip = 1, npop
         call report_column('ABS', m%channel_name(ip), pabs_pop(:,ip)*SCALE, &
                            t(1+ip,:), ok)
         call report_column('SCA', m%channel_name(ip), psca_pop(:,ip)*SCALE, &
                            t(1+ntype+ip,:), ok)
      end do

      ! The size integral must be the sum of the two halves over the
      ! populations; a mismatch would mean Cpol and Cpol_ext were built from
      ! different tables.
      dev = 0.0_wp
      do jw = 1, nwave
         if (Cpolext(jw) /= 0.0_wp) &
            dev = max(dev, abs((sum(pabs_pop(jw,:)) + sum(psca_pop(jw,:))) &
                               / Cpolext(jw) - 1.0_wp))
      end do
      write(*,'(a,es10.2)') '     ABS + SCA against the size integral   ', dev
      if (dev > 1.0e-12_wp .or. dev /= dev) ok = .false.

      call check_alignment_override(m, ok)

      if (.not. ok) write(*,'(a)') '   -> FAIL'
      deallocate(t, Cext, Cabs, Csca, Cpolext, pabs_pop, psca_pop)
   end subroutine compare_against_dustem


   subroutine check_alignment_override(m, ok)
      ! A host that installs its own alignment must override what the model
      ! file put there, on every aligned population, exactly as it does for
      ! astrodust. The three parameters below are not this model's and not the
      ! defaults, so an efficiency that still matches the DustEM law would show
      ! up as a mismatch here.
      type(dust_model_t), intent(inout) :: m
      logical,            intent(inout) :: ok
      real(wp), parameter :: FMAX = 0.8_wp, AAL = 0.05_wp, ALPHA = 2.5_wp
      real(wp) :: dev, f_want
      integer  :: ip, ja, st, naligned

      call dust_set_alignment(m, FMAX, AAL, ALPHA, status=st)
      ! status 4 only reports that no aligned scattering table matches, which
      ! is not loaded here; anything else is a rejection.
      if (st /= 0 .and. st /= 4) then
         write(*,'(a,i0)') '     dust_set_alignment refused, status = ', st
         ok = .false.;  return
      end if

      dev = 0.0_wp;  naligned = 0
      do ip = 1, size(m%pops)
         if (.not. allocated(m%pops(ip)%Cpol)) cycle
         naligned = naligned + 1
         do ja = 1, size(m%pops(ip)%aeff)
            f_want = falign_powerlaw(m%pops(ip)%aeff(ja), FMAX, AAL, ALPHA)
            dev = max(dev, abs(m%pops(ip)%falign(ja) - f_want))
         end do
      end do

      write(*,'(a,i0,a,es10.2)') '     dust_set_alignment refilled ', naligned, &
           ' aligned populations, max |f - power law| ', dev
      if (naligned == 0 .or. dev > 1.0e-14_wp .or. dev /= dev) ok = .false.
      if (m%align_tabulated) then
         write(*,'(a)') '     align_tabulated stayed .true. after dust_set_alignment'
         ok = .false.
      end if
   end subroutine check_alignment_override


   subroutine check_unaligned_model(label, grain_path, data_dir, ok)
      ! A model whose GRAIN file declares no POL population must come out with
      ! no polarized optics and a polarized extinction of exactly zero.
      character(len=*), intent(in)  :: label, grain_path, data_dir
      logical,          intent(out) :: ok
      type(dust_model_t)    :: m
      real(wp), allocatable :: Cext(:), Cabs(:), Csca(:), Cpolext(:)
      integer :: status

      ok = .true.
      write(*,'(a)') ' === '//label//', which declares no aligned population'
      call build_dustem(m, grain_path, data_dir, NT_IN, T_LO, T_HI, status, &
                        include_euv=.true.)
      if (status /= 0) then
         write(*,'(a,i0)') '   build_dustem failed, status = ', status
         ok = .false.;  return
      end if
      if (dust_has_polarized_optics(m)) then
         write(*,'(a)') '   it reports polarized optics'
         ok = .false.
      end if
      if (m%align_tabulated) then
         write(*,'(a)') '   align_tabulated is .true. for a model with no alignment law'
         ok = .false.
      end if
      allocate(Cext(m%NLAM), Cabs(m%NLAM), Csca(m%NLAM), Cpolext(m%NLAM))
      call size_integrated_extinction(m, Cext, Cabs, Csca, Cpol_ext=Cpolext)
      if (maxval(abs(Cpolext)) /= 0.0_wp) then
         write(*,'(a,es10.2)') '   polarized extinction is not zero, it reaches ', &
              maxval(abs(Cpolext))
         ok = .false.
      else
         write(*,'(a)') '   no polarized optics, polarized extinction exactly zero'
      end if
      deallocate(Cext, Cabs, Csca, Cpolext)
      if (.not. ok) write(*,'(a)') '   -> FAIL'
   end subroutine check_unaligned_model


   subroutine read_dustem_res(ref_path, ntype, nwave, t, ok)
      ! '#' header, a "ntype nwave" line, then nwave rows of 2*ntype+2 reals.
      character(len=*),      intent(in)  :: ref_path
      integer,               intent(out) :: ntype, nwave
      real(wp), allocatable, intent(out) :: t(:,:)
      logical,               intent(out) :: ok
      integer :: u, ios, i
      character(len=512) :: line

      ok = .true.;  ntype = 0;  nwave = 0
      open(newunit=u, file=ref_path, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(a)') '   cannot open '//ref_path;  ok = .false.;  return
      end if
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) then
            write(*,'(a)') '   no "ntype nwave" line in '//ref_path
            close(u);  ok = .false.;  return
         end if
         line = adjustl(line)
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#')   cycle
         read(line,*,iostat=ios) ntype, nwave
         if (ios /= 0) then
            write(*,'(a)') '   cannot parse "ntype nwave"'
            close(u);  ok = .false.;  return
         end if
         exit
      end do
      allocate(t(2*ntype+2, nwave))
      do i = 1, nwave
         read(u,*,iostat=ios) t(:, i)
         if (ios /= 0) then
            write(*,'(a,i0)') '   short reference file at row ', i
            close(u);  deallocate(t);  ok = .false.;  return
         end if
      end do
      close(u)
   end subroutine read_dustem_res


   subroutine report_column(what, name, mine, ref, ok)
      ! One column of the reference against ours.  A column the reference
      ! holds identically at zero is not a relative comparison and must not
      ! be made into one: the PAH population of G18 Model D is unaligned, so
      ! both its columns are zero, and dividing by them would give a NaN that
      ! then passes every tolerance test, NaN comparing false against
      ! anything.  Such a column is checked for being zero on our side too and
      ! is reported as zero rather than as an agreement.
      character(len=*), intent(in)    :: what, name
      real(wp),         intent(in)    :: mine(:), ref(:)
      logical,          intent(inout) :: ok
      real(wp) :: dev, scale_mine
      integer  :: jw, ncmp

      dev = 0.0_wp;  ncmp = 0
      do jw = 1, size(ref)
         if (ref(jw) /= 0.0_wp) then
            dev = max(dev, abs(mine(jw)/ref(jw) - 1.0_wp))
            ncmp = ncmp + 1
         end if
      end do

      if (ncmp == 0) then
         scale_mine = maxval(abs(mine))
         if (scale_mine == 0.0_wp) then
            write(*,'(a,a5,1x,a,a)') '     ', what, padded(name), &
                 '   zero in both, as this population is unaligned'
         else
            write(*,'(a,a5,1x,a,a,es10.2)') '     ', what, padded(name), &
                 '   reference is zero but ours reaches ', scale_mine
            ok = .false.
         end if
         return
      end if

      ! NaN fails: it compares false against every bound, so test it directly.
      if (dev /= dev) then
         write(*,'(a,a5,1x,a,a)') '     ', what, padded(name), &
              '        NaN in the comparison'
         ok = .false.
         return
      end if

      write(*,'(a,a5,1x,a,a,es10.2,a,i0,a)') '     ', what, padded(name), &
           '  ', dev, '   (', ncmp, ' wavelengths)'
      if (dev > TOL) ok = .false.
   end subroutine report_column


   pure function padded(name) result(s)
      ! The name right-justified in a fixed field, except that a name longer
      ! than the field widens it instead of losing characters.  An `a20` edit
      ! descriptor on a name declared longer than 20 prints the first twenty
      ! columns, which for a right-justified name are all blanks.
      character(len=*), intent(in) :: name
      integer, parameter :: NAMEW = 20
      character(len=max(NAMEW, len_trim(name))) :: s
      s = ''
      s(len(s)-len_trim(name)+1:) = trim(name)
   end function padded

end program test_dustem_polarized_extinction
