program test_dustem_extinction
   ! Compare build_dustem on each DustEM-defined model against DustEM's own
   ! extinction for the very same input files:
   !
   !   THEMIS        ../data/themis/reference/EXT_J13.RES        as distributed
   !   G18 Model D   ../data/g18d/reference/EXT_G17_ModelD.RES   see where.txt
   !
   ! Each reference is a DustEM run on the GRAIN file, optics tables and
   ! calorimetry this tree builds the model from, so the comparison isolates
   ! whether this tree reproduces DustEM's conventions -- the size-distribution
   ! formula, its mass normalization, the trapezoid in ln a the size sum has to
   ! be, and the linear-in-radius interpolation of the optics tables.  It is
   ! not a comparison of two dust models, and it never touches the solver:
   ! extinction is the size integral alone.
   !
   ! Reference format: '#' header, then a line "ntype nwave", then nwave rows of
   !   lambda[um]  ABS(1..ntype)  SCA(1..ntype)  EXT(total)
   ! with every cross section in cm^2 per H multiplied by 1e21.  The
   ! wavelengths are LAMBDA.DAT, so each model is built on the whole DustEM
   ! grid (include_euv) and no interpolation enters the comparison.
   !
   ! Passes when, for both models, the total extinction and every population's
   ! absorption and scattering agree to 1e-4 relative -- far inside the
   ! reference files' own 7-digit print precision, which they reach.
   use constants,         only: wp
   use sed_astrodust_mod, only: dust_model_t, build_dustem, size_integrated_extinction
   implicit none

   integer,  parameter :: NMODEL = 2
   integer,  parameter :: NT_IN  = 100
   real(wp), parameter :: T_LO = 1.0_wp, T_HI = 3000.0_wp
   real(wp), parameter :: TOL = 1.0e-4_wp
   real(wp), parameter :: SCALE = 1.0e21_wp

   character(len=16), parameter :: LABEL(NMODEL) = &
      [character(len=16) :: 'THEMIS', 'G18 Model D']
   character(len=64), parameter :: GRAIN(NMODEL) = &
      [character(len=64) :: '../data/themis/GRAIN_J13.DAT', &
                            '../data/g18d/GRAIN_G17_ModelD.DAT']
   character(len=64), parameter :: DDIR(NMODEL) = &
      [character(len=64) :: '../data/themis/', '../data/g18d/']
   character(len=64), parameter :: REF(NMODEL) = &
      [character(len=64) :: '../data/themis/reference/EXT_J13.RES', &
                            '../data/g18d/reference/EXT_G17_ModelD.RES']

   integer :: im
   logical :: ok, all_ok

   all_ok = .true.
   do im = 1, NMODEL
      call compare_against_dustem(trim(LABEL(im)), trim(GRAIN(im)), &
                                  trim(DDIR(im)), trim(REF(im)), ok)
      all_ok = all_ok .and. ok
      write(*,'(a)') ''
   end do

   if (.not. all_ok) then
      write(*,'(a,es10.2)') ' FAIL: tolerance is ', TOL;  stop 1
   end if
   write(*,'(a)') ' PASS'

contains

   subroutine compare_against_dustem(label, grain_path, data_dir, ref_path, ok)
      character(len=*), intent(in)  :: label, grain_path, data_dir, ref_path
      logical,          intent(out) :: ok

      type(dust_model_t)    :: m
      real(wp), allocatable :: t(:,:), Cext(:), Cabs(:), Csca(:)
      real(wp), allocatable :: abs_pop(:,:), sca_pop(:,:)
      real(wp) :: dev
      integer  :: u, ios, ntype, nwave, i, ip, ja, jw, status, npop
      character(len=512) :: line

      ok = .true.
      write(*,'(a)') ' === '//label

      ! ---- the DustEM reference -----------------------------------------
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

      ! ---- the model, on the whole DustEM grid ---------------------------
      call build_dustem(m, grain_path, data_dir, NT_IN, T_LO, T_HI, status, &
                        include_euv=.true.)
      if (status /= 0) then
         write(*,'(a,i0)') '   build_dustem failed, status = ', status
         deallocate(t);  ok = .false.;  return
      end if
      npop = size(m%pops)
      write(*,'(a,i0,a,i0,a)') '   ', npop, ' populations, ', m%NLAM, ' wavelengths'

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

      dev = 0.0_wp
      do jw = 1, nwave
         dev = max(dev, abs(m%lam(jw)/t(1,jw) - 1.0_wp))
      end do
      write(*,'(a,es10.2)') '   wavelength grid, max |model/reference - 1| : ', dev
      if (dev > 1.0e-6_wp) ok = .false.

      ! ---- the size integrals --------------------------------------------
      allocate(Cext(nwave), Cabs(nwave), Csca(nwave))
      call size_integrated_extinction(m, Cext, Cabs, Csca)

      allocate(abs_pop(nwave, npop), sca_pop(nwave, npop))
      abs_pop = 0.0_wp;  sca_pop = 0.0_wp
      do ip = 1, npop
         do ja = 1, size(m%pops(ip)%dn)
            do jw = 1, nwave
               abs_pop(jw, ip) = abs_pop(jw, ip) + m%pops(ip)%dn(ja) * m%pops(ip)%Cabs(jw, ja)
               sca_pop(jw, ip) = sca_pop(jw, ip) + m%pops(ip)%dn(ja) * m%pops(ip)%Csca(jw, ja)
            end do
         end do
      end do

      ! ---- compare --------------------------------------------------------
      write(*,'(a)') '   max |model/reference - 1|'
      call report_column('total extinction', '', Cext*SCALE, t(2*ntype+2,:), ok)
      do ip = 1, npop
         call report_column('ABS', m%channel_name(ip), abs_pop(:,ip)*SCALE, &
                            t(1+ip,:), ok)
         call report_column('SCA', m%channel_name(ip), sca_pop(:,ip)*SCALE, &
                            t(1+ntype+ip,:), ok)
      end do
      if (.not. ok) write(*,'(a)') '   -> FAIL'
      deallocate(t, Cext, Cabs, Csca, abs_pop, sca_pop)
   end subroutine compare_against_dustem


   subroutine report_column(what, name, mine, ref, ok)
      ! One column of the reference against ours.  A column the reference
      ! holds identically at zero is not a relative comparison and must not
      ! be made into one: DustEM's PAH table for G18 Model D carries no
      ! scattering at all, and dividing by it produced a NaN that then passed
      ! every tolerance test, NaN comparing false against anything.  Such a
      ! column is checked for being zero on our side too and is reported as
      ! zero rather than as an agreement.
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
         ! The reference is zero everywhere; ours has to be as well.
         scale_mine = maxval(abs(mine))
         if (scale_mine == 0.0_wp) then
            write(*,'(a,a5,1x,a,a)') '     ', what, padded(name), &
                 '   zero in both, as the table has no scattering'
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

end program test_dustem_extinction
