module dust_lib
   ! RT-facing entry point for the model-agnostic dust thermal-emission
   ! library. A Fortran 3D radiative-transfer code links this module and:
   !
   !   use dust_lib
   !   type(dust_model_t) :: m
   !   call build_astrodust(m, qtab, sizedist, NT, T_lo, T_hi)   ! once
   !   ...
   !   do icell = 1, ncells
   !      ! ... assemble local mean intensity J_lam(:) on m%lam ...
   !      call dust_emission(m, J_lam, lamI_total [, lamI_chan])  ! per cell
   !   end do
   !
   ! EXTREME ULTRAVIOLET. The astrodust model's wavelength grid is the T-matrix
   ! Q table's, which stops at 0.0912 um (13.6 eV, the Lyman limit). A host
   ! that transports harder photons passes the shortest wavelength it needs as
   ! the optional lam_min [um]:
   !   call build_astrodust(m, qtab, sizedist, NT, T_lo, T_hi, lam_min=0.0124_wp)
   ! The grid is then carried from lam_min up to the table, log-spaced no more
   ! coarsely than the table's own dln(lambda) = 0.01156, and m%lam is that
   ! longer grid; everything else about the API is unchanged. OMITTING lam_min
   ! leaves the grid, and every number the model produces, exactly as before.
   !
   ! In the extended band the astrodust optics are computed from the DH21
   ! dielectric function by the T-MATRIX, on the same b/a = 1.400 oblate
   ! spheroid and the same random-orientation average as the Q table itself
   ! (spheroid_q, which takes the Rayleigh dipole limit below x = 0.1 and the
   ! geometric-optics limit above x = 50, exactly as the table's own driver
   ! does). The band is therefore a continuation of the table, not a different
   ! particle glued to it: measured at the seam, C_abs matches a log-log
   ! extrapolation of the table to 0.07% at worst over the whole size range,
   ! and <cos> agrees exactly where the table's own regime bounds put it at 0.
   !
   ! That spheroid is the ONLY thing in this library that needs the T-matrix,
   ! and it is not compiled in by default. The plain libsedust.a therefore has
   ! no T-matrix in it and links without one; asking it for the spheroid
   ! (euv_tmatrix = .true., the default whenever lam_min is given) is REFUSED
   ! with status 6 rather than answered with a sphere. To get the spheroid:
   !
   !   cd tmatrix && make libtmatrix.a
   !   cd sed     && WITH_TMATRIX=1 make libsedust.a     # or ./build_lib.sh
   !
   ! and, once in the host, before the model is built:
   !
   !   use euv_astrodust_tmatrix, only: use_tmatrix_euv_band_optics
   !   call use_tmatrix_euv_band_optics()
   !
   ! (Or register an implementation of your own: sed_register_euv_band_optics
   ! takes any procedure matching the abstract interface euv_band_optics_i,
   ! both re-exported here, and sed_forget_euv_band_optics undoes it.)
   !
   ! That costs real time. Measured here on 167 radii (gfortran 13.1, -O3),
   ! band only:
   !
   !   lam_min                 n_euv  band pts  T-matrix pts  threads    time
   !   0.0247968 um (50 eV)      113    18,871        12,204        1   527.5 s
   !                                                               2   264.2 s
   !                                                               4   136.7 s
   !                                                               8    69.2 s
   !   1.001e-4 um (12.4 keV)    590    98,530        41,900       72    73.5 s
   !
   ! so the band is built with OpenMP (one T-matrix workspace for each thread,
   ! 80 MiB apiece once its full-direct arrays are allocated; nothing is
   ! allocated when there is no EUV band). Peak resident memory grows with it:
   ! 64 / 110 / 197 / 376 MB at 1 / 2 / 4 / 8 threads, 3.0 GB at 72. The size
   ! loop carries no dependence, so the result does not depend on the thread
   ! count or the schedule -- the four runs above were byte-identical.
   !
   ! A host that cannot afford this passes euv_tmatrix = .false.:
   !   call build_astrodust(m, qtab, sizedist, NT, T_lo, T_hi, &
   !                        lam_min=0.0124_wp, euv_tmatrix=.false.)
   ! which substitutes Bohren-Huffman Mie for the volume-equivalent SPHERE and
   ! runs in a tenth of a second. It is an APPROXIMATION: good to ~2% where the
   ! T-matrix regime holds, but converging on a fixed ~2% underestimate for
   ! grains in the geometric-optics limit rather than on zero, and leaving a
   ! visible step at the seam (median 1.2%, worst 16% in C_abs; <cos> jumps
   ! from 0.9 to 0 for a >~ 1 um). q_astrodust.f90 documents it size by size.
   !
   ! The carbonaceous optics have no such seam: above the DL07 PAH cutoff
   ! (21.4 eV) they are Mie on the D03 graphite dielectric functions, which
   ! run to 6.2e-5 um.
   !
   ! The dielectric function must be the one the Q table was computed from --
   ! same porosity, iron fraction and axial ratio -- or the model changes
   ! material at the seam. build_astrodust takes it as the optional
   ! astrodust_index_path for that reason, and reports the file whenever the
   ! extension is active; omitted, it is the P = 0.20, f_Fe = 0.00, b/a = 1.400
   ! file that pairs with q_astrodust_P0.20_Fe0.00_1.400.dat. lam_min shorter
   ! than that file's own shortest wavelength (1.00003e-4 um = 12.4 keV) is
   ! refused rather than served with a frozen refractive index.
   !
   ! Two usage modes:
   !   (a) single-cell EXACT solve: dust_emission(m, J_lam, ...)  -- arbitrary J(lambda).
   !   (b) precomputed TABLE + interpolation, for when the cell field is a
   !       fixed reference SHAPE scaled by an intensity U:
   !         call dust_build_table(m, J_ref, U_grid, tab)        ! once
   !         call dust_emission_interp(tab, U, lamI_total [, lamI_chan]) ! per cell
   !
   ! NONLINEARITY CAVEAT: stochastic heating makes the emission a NONLINEAR
   ! functional of the full J(lambda), not of a scalar. The table is valid
   ! ONLY when the cell field is  U * (fixed J_ref shape).  For cells whose
   ! field SHAPE departs from J_ref (e.g. hardened spectra near hot stars),
   ! use the single-cell exact dust_emission instead. The table is built from
   ! exact solves, so it reproduces the U grid points to round-off, except that
   ! dust_emission_interp floors the stored emissivities at 1e-300 before the
   ! log, so a node whose true value is 0 comes back as 1e-300, not 0;
   ! interpolation between grid points is the usual smooth-in-U approximation.
   !
   ! dust_emission takes an optional final argument, status (integer):
   !   call dust_emission(m, J_lam, lamI_total [, lamI_chan] [, status])
   ! On return status = 0 means success, 1 an unknown m%stoch_method, and 2 a
   ! 'qm' model whose populations are missing their radii. When status is
   ! present a bad model is reported through it instead of stopping the
   ! process; when it is omitted such a model stops the run, as before.
   !
   ! EXTINCTION. dust_extinction is the extinction counterpart of
   ! dust_emission, so an RT host takes its opacity from the same model object
   ! and on the same wavelength grid (m%lam) as its emission:
   !   call dust_extinction(m, Cext, Cabs, Csca [, gbar] [, albedo] [, status])
   ! All three required outputs are (m%NLAM) cross sections per H atom
   ! [cm^2/H], with Cext = Cabs + Csca. Optional gbar is the
   ! scattering-weighted asymmetry <cos>; optional albedo is Csca/Cext, 0 where
   ! Cext underflows.
   !
   ! What it returns is the PRECOMPUTED size-integrated curve of a
   ! data/kext_*.dat table, interpolated onto m%lam -- log(lambda)-log(C) for
   ! the cross sections, linear in log(lambda) for <cos>, and the table value
   ! itself wherever a model wavelength coincides with a table wavelength.
   ! Nothing is integrated over grain size at call time. The size integral is
   ! size_integrated_extinction, which is what wrote those tables and is what
   ! the standalone calc_kext.x calls; it is re-exported here for a host that
   ! wants to redo the integral rather than read the product.
   !
   ! The BUILDER loads the table, not dust_extinction, because the table paths
   ! are relative to the sed/ directory and a host is free to change directory
   ! once the model is built. Each builder takes an optional kext_path:
   !   call build_astrodust(m, qtab, sizedist, NT, T_lo, T_hi, &
   !                        kext_path='../data/kext_astrodust_MW.dat')
   ! A file named that way is a contract: if it cannot be read the build FAILS.
   ! Omitting it falls to the model's default table --
   !   astrodust   ../data/kext_astrodust_MW_euv.dat
   !   dl07        ../data/kext_dl07_MW_euv.dat
   !   zubko       ../data/kext_zubko_BARE_GR_S.dat
   !   from_files  none
   ! -- and a default that cannot be read leaves the model with no extinction
   ! to serve (dust_extinction then reports status 2) but does NOT fail the
   ! build, so a driver that only wants emission still runs.
   !
   ! The defaults are the EUV products because their wavelength range CONTAINS
   ! the plain grid's: the two are the same size integral on the same optics,
   ! the extended one merely carried below the T-matrix Q table's 0.0912 um
   ! end, and their nodes coincide over the overlap. One default therefore
   ! covers a host that transports ionizing radiation and one that does not.
   ! It is the NARROWER non-EUV table that has to be named explicitly. A model
   ! grid reaching outside the table it was given is refused (status 3), never
   ! extrapolated.
   !
   ! Every model carries its scattering optics, so all four write a physical
   ! albedo into their tables:
   !   astrodust  Csca and <cos> from the random-orientation T-matrix Q table
   !   dl07       Mie on the D03 astrosilicate and graphite dielectric
   !              functions (q_silicate_full / q_graphite_full); the PAH
   !              component scatters negligibly and enters through absorption
   !   zubko      Q_sca and <cos> columns of the model's own DustEM tables
   !   from_files same, for whatever tables the descriptor names
   ! A population that genuinely does not scatter leaves its optics unallocated
   ! and contributes zero to the size integral. dust_extinction's status is 0
   ! on success, 1 if an output array is not of size m%NLAM, 2 if no table was
   ! loaded, and 3 if m%lam reaches outside the table; when status is omitted
   ! such a call stops the run.
   !
   ! DUST MASS. The cross sections above are per H nucleon. A host that carries
   ! a dust mass density instead converts with the model's own dust mass per H:
   !   Mdust_H = dust_mass_per_H(m)                ! [g/H]
   !   K_abs   = Cabs / Mdust_H                    ! [cm^2/g]
   ! It is one wavelength-independent constant, summed from the model's size
   ! distribution and the solid density of each population's material,
   ! M_dust/N_H = sum_pop rho_bulk * sum_a (4/3) pi a^3 dn(a). Every builder
   ! states those densities, so it is defined for every model in the tree:
   !   astrodust  rho_Ad = 2.74, rho_PAH = 2.0 g/cm^3 (HD23)
   !   dl07       astrosilicate 3.5, graphite 2.24 g/cm^3 (Draine & Lee 1984,
   !              Weingartner & Draine 2001)
   !   zubko      each component's density as its own DustEM optics file states it
   !   from_files the descriptor's rho, or the optics file's when that is 0
   ! It is the same number calc_kext.x normalizes the trailing K_abs column of
   ! every data/kext_*.dat by.
   !
   ! dust_build_table and dust_emission_interp take the same optional final
   ! status argument (0 = success); when present a bad argument is reported
   ! through it instead of stopping the process; when omitted such a call stops
   ! the run.
   !   dust_build_table:      1 size(J_ref) /= m%NLAM
   !                          2 size(U_grid) < 2
   !                          3 U_grid not positive-and-strictly-increasing
   !   dust_emission_interp:  1 U <= 0
   !                          2 size(lamI_total) /= tab%NLAM
   !                          3 lamI_chan present but not (tab%NLAM, tab%n_channel)
   !
   ! The model builders take the same optional final argument, status (integer,
   ! 0 = success). When present, a missing or malformed input file is reported
   ! through it and the model is NOT built (so an RT host can recover); when
   ! omitted such a failure stops the process, which the CLI drivers rely on.
   ! The exact non-zero code only distinguishes the failing stage; the contract
   ! is simply "0 = built, non-zero = build failed". Codes per builder:
   !   build_astrodust / build_dl07:  1 Q-table load failed
   !                                  2 size-distribution load failed
   !                                  5 an explicitly named kext_path failed to
   !                                    load
   !   build_astrodust only, and only when lam_min asks for the EUV band:
   !                                  3 astrodust dielectric function load failed
   !                                  4 lam_min below that dielectric function's
   !                                    shortest wavelength
   !                                  6 euv_tmatrix = .true. but the spheroid
   !                                    optics of that band are not available:
   !                                    no euv_band_optics_i is registered (a
   !                                    library built without the T-matrix), or
   !                                    the registered one could not compute
   !                                    the band
   !   build_zubko:   1 config read failed        2 fewer than 3 components
   !                  3 a component's optics read  4 grid inconsistent
   !                  5 a component's calorimetry read failed
   !                  6 an explicitly named kext_path failed to load
   !   build_from_files: 1 descriptor open   2 too many pop: lines
   !                     3 invalid channel   4 no pop: lines
   !                     5 optics read       6 grid inconsistent
   !                     7 size-dist read    8 calorimetry read failed
   !                     9 an explicitly named kext_path failed to load
   !
   ! LINKING. One archive, and one .mod search path:
   !   cd sed && make libsedust.a          # or ./build_lib.sh -> sed/lib/
   !   gfortran -I<sed> my_rt.f90 -L<sed> -lsedust -fopenmp -o my_rt.x
   ! Nothing in this library reaches the T-matrix, so ../tmatrix need not even
   ! be built.
   !
   ! A host that wants the astrodust EUV band on the spheroid builds the
   ! T-matrix into the SAME archive, so the link line does not change:
   !   cd tmatrix && make libtmatrix.a
   !   cd sed     && WITH_TMATRIX=1 make libsedust.a
   !   gfortran -I<sed> my_rt.f90 -L<sed> -lsedust -fopenmp -o my_rt.x
   ! and calls use_tmatrix_euv_band_optics() once, as above. libtmatrix.a
   ! itself needs no OpenMP runtime; the threading is on the SEDust side.
   !
   ! The validated solver core (sed_grain_loop & helpers in sed_astrodust_mod)
   ! is untouched; this module only re-exports the model API and adds the
   ! table/interpolation layer.
   use constants,         only: wp
   use sed_mathlib,           only: locate
   use sed_astrodust_mod, only: dust_model_t, &
                                build_astrodust, build_dl07, build_zubko, build_from_files, &
                                dust_emission, dust_emission_single_teq, &
                                dust_extinction, size_integrated_extinction, &
                                dust_mass_per_H, euv_band_optics_i, &
                                sed_register_euv_band_optics, &
                                sed_forget_euv_band_optics
   implicit none
   private

   ! Re-exported model API
   public :: dust_model_t, build_astrodust, build_dl07, build_zubko, build_from_files
   public :: dust_emission, dust_emission_single_teq, dust_extinction
   public :: size_integrated_extinction, dust_mass_per_H
   ! Optics of the astrodust EUV band; see the EXTREME ULTRAVIOLET note above.
   public :: euv_band_optics_i
   public :: sed_register_euv_band_optics, sed_forget_euv_band_optics
   ! Table API
   public :: dust_emis_table_t, dust_build_table, dust_emission_interp, dust_free_table
   ! Convenience accessors
   public :: dust_nlam, dust_lambda, dust_n_channel, dust_channel_name

   type :: dust_emis_table_t
      integer               :: NLAM = 0, n_channel = 0, NU = 0
      real(wp), allocatable :: U(:)            ! (NU)   intensity-scaling grid
      real(wp), allocatable :: logU(:)         ! (NU)   log(U), cached for interp
      real(wp), allocatable :: lam(:)          ! (NLAM) [um]
      real(wp), allocatable :: J_ref(:)        ! (NLAM) reference field shape (U=1)
      real(wp), allocatable :: total(:,:)      ! (NLAM, NU)
      real(wp), allocatable :: chan(:,:,:)     ! (NLAM, n_channel, NU)
   end type dust_emis_table_t

contains

   ! --- accessors -------------------------------------------------------
   pure integer function dust_nlam(m)
      type(dust_model_t), intent(in) :: m
      dust_nlam = m%NLAM
   end function dust_nlam

   pure integer function dust_n_channel(m)
      type(dust_model_t), intent(in) :: m
      dust_n_channel = m%n_channel
   end function dust_n_channel

   function dust_lambda(m) result(lam)
      type(dust_model_t), intent(in) :: m
      real(wp), allocatable :: lam(:)
      lam = m%lam
   end function dust_lambda

   function dust_channel_name(m, ic) result(name)
      type(dust_model_t), intent(in) :: m
      integer,            intent(in) :: ic
      character(len=16) :: name
      name = m%channel_name(ic)
   end function dust_channel_name

   ! --- emission table over an intensity grid ---------------------------
   subroutine dust_build_table(m, J_ref, U_grid, tab, status)
      ! Precompute lamI(lambda) for J = U*J_ref at each U in U_grid (which must
      ! be positive and strictly ascending). m must be the active model.
      type(dust_model_t),      intent(in)  :: m
      real(wp),                intent(in)  :: J_ref(:)    ! (NLAM) shape at U=1
      real(wp),                intent(in)  :: U_grid(:)   ! (NU)
      type(dust_emis_table_t), intent(out) :: tab
      ! Optional error report (0 = success); see the module header for codes.
      integer, optional,       intent(out) :: status
      real(wp), allocatable :: total(:), chan(:,:)
      integer :: iu, nu

      if (present(status)) status = 0
      nu = size(U_grid)

      if (size(J_ref) /= m%NLAM) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,i0,a,i0)') 'dust_build_table: size(J_ref)=', size(J_ref), &
                                    ' /= m%NLAM=', m%NLAM
            stop 1
         end if
      end if
      if (nu < 2) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a,i0)') 'dust_build_table: need size(U_grid) >= 2, got ', nu
            stop 1
         end if
      end if
      if (U_grid(1) <= 0.0_wp .or. any(U_grid(2:nu) <= U_grid(1:nu-1))) then
         if (present(status)) then
            status = 3;  return
         else
            write(*,'(a)') 'dust_build_table: U_grid must be positive and strictly increasing'
            stop 1
         end if
      end if

      call dust_free_table(tab)
      tab%NLAM = m%NLAM;  tab%n_channel = m%n_channel;  tab%NU = size(U_grid)
      allocate(tab%U(tab%NU), tab%logU(tab%NU), tab%lam(tab%NLAM), tab%J_ref(tab%NLAM))
      allocate(tab%total(tab%NLAM, tab%NU))
      allocate(tab%chan(tab%NLAM, tab%n_channel, tab%NU))
      tab%U = U_grid;  tab%lam = m%lam;  tab%J_ref = J_ref
      tab%logU = log(tab%U)   ! cached once; the U bracket is grid-fixed at interp

      allocate(total(m%NLAM), chan(m%NLAM, m%n_channel))
      do iu = 1, tab%NU
         call dust_emission(m, U_grid(iu)*J_ref, total, chan)
         tab%total(:, iu)  = total
         tab%chan(:, :, iu) = chan
      end do
      deallocate(total, chan)
   end subroutine dust_build_table

   subroutine dust_emission_interp(tab, U, lamI_total, lamI_chan, status)
      ! Log-log interpolate the table at intensity U (per wavelength and
      ! channel). Reproduces the U grid points to round-off, except that the
      ! stored emissivities are floored at 1e-300 before the log, so a node
      ! whose true value is 0 comes back as 1e-300, not 0; clamps to the grid
      ! ends outside [U(1), U(NU)].
      type(dust_emis_table_t), intent(in)  :: tab
      real(wp),                intent(in)  :: U
      real(wp),                intent(out) :: lamI_total(:)      ! (NLAM)
      real(wp), optional,      intent(out) :: lamI_chan(:,:)     ! (NLAM, n_channel)
      ! Optional error report (0 = success); see the module header for codes.
      ! Validation is scalar/size-only to keep this on the hot path.
      integer, optional,       intent(out) :: status
      real(wp) :: ly(tab%NU)
      real(wp) :: lr, lUq
      integer  :: k, c, jlo

      if (present(status)) status = 0
      if (U <= 0.0_wp) then
         if (present(status)) then
            status = 1;  return
         else
            write(*,'(a,es12.4)') 'dust_emission_interp: need U > 0, got ', U
            stop 1
         end if
      end if
      if (size(lamI_total) /= tab%NLAM) then
         if (present(status)) then
            status = 2;  return
         else
            write(*,'(a,i0,a,i0)') 'dust_emission_interp: size(lamI_total)=', &
                                    size(lamI_total), ' /= tab%NLAM=', tab%NLAM
            stop 1
         end if
      end if
      if (present(lamI_chan)) then
         if (size(lamI_chan,1) /= tab%NLAM .or. size(lamI_chan,2) /= tab%n_channel) then
            if (present(status)) then
               status = 3;  return
            else
               write(*,'(a)') 'dust_emission_interp: lamI_chan must be (tab%NLAM, tab%n_channel)'
               stop 1
            end if
         end if
      end if

      ! tab%logU (= log(tab%U)) is cached at build time and is strictly
      ! ascending because U_grid is positive and strictly increasing. The U
      ! bracket is the same for every wavelength and channel, so locate it once
      ! and reuse jlo below. The interpolation arithmetic is the ascending
      ! branch of interp1, kept term-for-term (log-log, clamped to the ends).
      lUq = log(U)
      call locate(tab%logU, lUq, jlo)
      do k = 1, tab%NLAM
         ly = log(max(tab%total(k, :), 1.0e-300_wp))
         if (jlo == 0) then
            lr = ly(1)
         else if (jlo == tab%NU) then
            lr = ly(tab%NU)
         else
            lr = ly(jlo) + (ly(jlo+1)-ly(jlo))*(lUq-tab%logU(jlo))/(tab%logU(jlo+1)-tab%logU(jlo))
         end if
         lamI_total(k) = exp(lr)
      end do
      if (present(lamI_chan)) then
         do c = 1, tab%n_channel
            do k = 1, tab%NLAM
               ly = log(max(tab%chan(k, c, :), 1.0e-300_wp))
               if (jlo == 0) then
                  lr = ly(1)
               else if (jlo == tab%NU) then
                  lr = ly(tab%NU)
               else
                  lr = ly(jlo) + (ly(jlo+1)-ly(jlo))*(lUq-tab%logU(jlo))/(tab%logU(jlo+1)-tab%logU(jlo))
               end if
               lamI_chan(k, c) = exp(lr)
            end do
         end do
      end if
   end subroutine dust_emission_interp

   subroutine dust_free_table(tab)
      type(dust_emis_table_t), intent(inout) :: tab
      if (allocated(tab%U))     deallocate(tab%U)
      if (allocated(tab%logU))  deallocate(tab%logU)
      if (allocated(tab%lam))   deallocate(tab%lam)
      if (allocated(tab%J_ref)) deallocate(tab%J_ref)
      if (allocated(tab%total)) deallocate(tab%total)
      if (allocated(tab%chan))  deallocate(tab%chan)
      tab%NLAM = 0;  tab%n_channel = 0;  tab%NU = 0
   end subroutine dust_free_table

end module dust_lib
