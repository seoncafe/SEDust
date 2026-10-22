! Fifteen fixed radiation fields from the original timing driver. The shape
! normalization always uses the original full-grid logarithmic trapezoid.
! Usage: validation.x 200 heuristic[,heuristic#1e-6] output_directory data_directory
program validation_driver
   use dust_lib
   use omp_lib, only: omp_get_wtime
#ifdef BIN_POWER
   use sed_astrodust_mod, only: emission_bin_power_cutoff, &
      emission_bin_power_total, emission_bin_power_skipped, &
      emission_bin_power_weighted_total, emission_bin_power_weighted_skipped
#endif
   implicit none
   integer, parameter :: dp=kind(1d0), nrep=5, nshape=3, npow=5
   real(dp), parameter :: pow(npow)=[3d-23,1d-22,1d-21,1d-20,1d-19]
   character(len=2), parameter :: shn(nshape)=['S0','S1','S2']
   real(dp), parameter :: hpl=6.62607015d-34, cl=2.99792458d8, kb=1.380649d-23
   real(dp), parameter :: pi=3.141592653589793d0, hc_evum=1.23984d0, tbb=4.0d4
   character(len=512) :: data_dir
   type(dust_model_t) :: m
   integer :: nt, st, nlam, nch, ish, ipw, irep, im, nmeth, u, k, cut
   character(len=256) :: arg, outdir, methods, meth(4), msg, tag, path, field_tag
   real(dp), allocatable :: lam(:), cext(:), cabs(:), csca(:), wln(:)
   real(dp), allocatable :: shape(:,:), jfield(:), lami(:), lami_ch(:,:)
   real(dp), allocatable :: prev(:), prev_ch(:,:)
   real(dp) :: p0(nshape), times(nrep), e, tau, scale, p_old, p_supported, pemit
   real(dp) :: raw_total, raw_skipped, weighted_total, weighted_skipped
   real(dp) :: t0, build_time
   logical :: identical

   call get_command_argument(1,arg); read(arg,*) nt
   call get_command_argument(2,methods)
   call get_command_argument(3,outdir)
   call get_command_argument(4,data_dir)
   if (len_trim(data_dir)==0) error stop 'data_directory argument is required'
   call split_csv(methods,meth,nmeth)
   if (nmeth<1 .or. nmeth>size(meth)) error stop 'invalid method count'
   t0=omp_get_wtime()
   call build_dust(m,'astrodust',data_dir,nt,2.7d0,5.0d3,include_euv=.true., &
      status=st,message=msg,sd_index=7,u_isrf=1.0d0)
   build_time=omp_get_wtime()-t0
   if (st/=0) then
      write(*,*) 'build_dust status=',st,trim(msg)
      stop 2
   end if
   nlam=dust_nlam(m); nch=dust_n_channel(m)
   lam=dust_lambda(m)
   allocate(cext(nlam),cabs(nlam),csca(nlam),wln(nlam))
   call dust_extinction(m,cext,cabs,csca,status=st)
   if (st/=0) stop 3
   wln(1)=0.5d0*log(lam(2)/lam(1))
   wln(nlam)=0.5d0*log(lam(nlam)/lam(nlam-1))
   do k=2,nlam-1
      wln(k)=0.5d0*log(lam(k+1)/lam(k-1))
   end do
   allocate(shape(nlam,nshape),jfield(nlam),lami(nlam),lami_ch(nlam,nch), &
      prev(nlam),prev_ch(nlam,nch))
   shape=0d0
   do k=1,nlam
      e=hc_evum/lam(k)
      if (e>=4.78d0 .and. e<=100d0) then
         shape(k,1)=planck_si(lam(k))
         tau=0d0
         if (e>=13.6d0) tau=(e/13.6d0)**(-3)
         shape(k,2)=planck_si(lam(k))*exp(-tau)
      end if
      if (e>=4.78d0 .and. e<13.6d0) shape(k,3)=planck_si(lam(k))
   end do
   do ish=1,nshape
      p0(ish)=4d0*pi*sum(cabs*shape(:,ish)*1d-3*lam*wln)
      if (p0(ish)<=0d0) stop 4
   end do
   open(newunit=u,file=trim(outdir)//'/lambda.bin',access='stream',status='replace')
   write(u) nlam,nch,lam,cabs
   close(u)
   open(newunit=u,file=trim(outdir)//'/timing.tsv',status='replace')
   write(u,'(a)') '# method shape power_target pabs_old pabs_supported pemit '// &
      'time_median time_min time_max identical raw_total raw_skipped '// &
      'weighted_total weighted_skipped'
   do im=1,nmeth
      cut=index(meth(im),'#')
      if (cut>0) then
         m%stoch_method=meth(im)(1:cut-1)
#ifdef BIN_POWER
         read(meth(im)(cut+1:),*) emission_bin_power_cutoff
#else
         error stop 'emission bin cutoff requires BIN_POWER build'
#endif
      else
         m%stoch_method=trim(meth(im))
#ifdef BIN_POWER
         emission_bin_power_cutoff=0d0
#endif
      end if
      do ish=1,nshape
         do ipw=1,npow
            scale=pow(ipw)/p0(ish)
            jfield=scale*shape(:,ish)
            write(field_tag,'(a,"_p",i0)') shn(ish),ipw
            path=trim(outdir)//'/field_'//trim(field_tag)//'.bin'
            open(newunit=k,file=trim(path),access='stream',status='replace')
            write(k) jfield
            close(k)
            p_old=4d0*pi*sum(cabs*jfield*1d-3*lam*wln)
#ifdef FIELD_SUPPORT
            p_supported=4d0*pi*1d-3*radiation_absorption(lam,cabs,jfield)
#else
            p_supported=supported_power(lam,cabs,jfield)
#endif
            identical=.true.
            raw_total=0d0; raw_skipped=0d0
            weighted_total=0d0; weighted_skipped=0d0
            do irep=1,nrep
               t0=omp_get_wtime()
               call dust_emission(m,jfield,lami,lami_ch,status=st)
               times(irep)=omp_get_wtime()-t0
               if (st/=0) then
                  write(*,*) 'dust_emission status=',st,' method=',trim(meth(im)), &
                     ' shape=',shn(ish),' power=',ipw
                  stop 5
               end if
               if (irep>1) identical=identical .and. all(lami==prev) .and. all(lami_ch==prev_ch)
               prev=lami; prev_ch=lami_ch
#ifdef BIN_POWER
               raw_total=emission_bin_power_total
               raw_skipped=emission_bin_power_skipped
               weighted_total=emission_bin_power_weighted_total
               weighted_skipped=emission_bin_power_weighted_skipped
#endif
               write(tag,'(a,"_NT",i0,"_",a,"_p",i0,"_r",i0)') &
                  trim(meth(im)),nt,shn(ish),ipw,irep
               path=trim(outdir)//'/spec_'//trim(tag)//'.bin'
               open(newunit=k,file=trim(path),access='stream',status='replace')
               write(k) lami,lami_ch
               close(k)
            end do
            pemit=4d0*pi*sum(lami*wln)
            write(u,'(a,1x,a,1x,7(es24.16,1x),l1,1x,4(es24.16,1x))') &
               trim(meth(im)),shn(ish),pow(ipw),p_old,p_supported,pemit, &
               median(times),minval(times),maxval(times),identical, &
               raw_total,raw_skipped,weighted_total,weighted_skipped
            flush(u)
            write(*,'(a,1x,a,1x,es12.4)') trim(meth(im)),trim(field_tag),median(times)
         end do
      end do
   end do
   close(u)
   write(*,'(a,i0,a,i0,a,es12.4)') 'VALIDATION_DONE NLAM=',nlam, &
      ' NCHANNEL=',nch,' BUILD_SECONDS=',build_time
contains
   function planck_si(lum) result(b)
      real(dp),intent(in) :: lum
      real(dp) :: b,lm,x
      lm=lum*1d-6
      x=hpl*cl/(lm*kb*tbb)
      b=2d0*hpl*cl**2/lm**5/(exp(x)-1d0)
   end function planck_si
   function supported_power(lambda, cross, field) result(power)
      real(dp),intent(in) :: lambda(:),cross(:),field(:)
      real(dp) :: power, floor_value, dln
      integer :: i,lo,hi,n
      n=size(lambda); power=0d0
      floor_value=1d-12*maxval(field)
      if (floor_value<=0d0) return
      lo=0; hi=0
      do i=1,n
         if (field(i)>floor_value .and. lambda(i)>0d0) then
            if (lo==0) lo=i
            hi=i
         end if
      end do
      if (lo==hi) return
      do i=lo,hi-1
         dln=log(lambda(i+1)/lambda(i))
         power=power+0.5d0*(cross(i)*field(i)*lambda(i)+ &
            cross(i+1)*field(i+1)*lambda(i+1))*dln
      end do
      power=4d0*pi*1d-3*power
   end function supported_power
   function median(a) result(md)
      real(dp),intent(in) :: a(:)
      real(dp) :: md,b(size(a)),t
      integer :: i,j,n
      b=a; n=size(b)
      do i=2,n
         t=b(i); j=i-1
         do while (j>=1)
            if (b(j)<=t) exit
            b(j+1)=b(j); j=j-1
         end do
         b(j+1)=t
      end do
      md=b((n+1)/2)
   end function median
   subroutine split_csv(s,out,n)
      character(len=*),intent(in) :: s
      character(len=*),intent(out) :: out(:)
      integer,intent(out) :: n
      integer :: start,i
      n=0; start=1
      do i=1,len_trim(s)+1
         if (i>len_trim(s)) then
            n=n+1
            if (n>size(out)) error stop 'too many emission methods'
            out(n)=s(start:i-1)
         else if (s(i:i)==',') then
            n=n+1
            if (n>size(out)) error stop 'too many emission methods'
            out(n)=s(start:i-1); start=i+1
         end if
      end do
   end subroutine split_csv
end program validation_driver
