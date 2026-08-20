module sed_apply_options
   ! Where a run setting LANDS: the module state the solver, the radiation
   ! field and the optics read.  It is kept apart from sed_run_options -- which
   ! owns the vocabulary, the parsing, the refusals and the filename tag --
   ! because those are the same for every program here, while this is not: a
   ! program that computes an enthalpy table or a cross-section table links no
   ! solver and no radiation field, and would otherwise have to drag both in to
   ! get the shared command line.
   use constants,              only: wp
   use radfield,               only: use_mathis_corrected, bbody
   use sed_astrodust_mod,      only: stoch_method, use_induced_emission, &
                                     gd_photon_cutoff
   use stoch_qm_mod,           only: qm_method, qm_nstate_default, qm_nisrf_max
   use qpah,                   only: qpah_graphite_source, qpah_xsec_vintage
   use enthalpy_astrodust_mod, only: s1_density_corrected
   use sed_run_options,        only: run_options_t, T_HARD_EUV, W_HARD_EUV, &
                                     LAM_LYMAN_UM
   implicit none
   private
   public :: apply_run_options, add_hard_euv_component

contains

   subroutine apply_run_options(o)
      ! Push the settings into the module state the solver reads.  Call it
      ! BEFORE the model is built: the radiation-field convention fixes the CMB
      ! temperature that goes into the cooling term, and the graphite source
      ! fixes the optics of the xi blend, both of which are built once.  A
      ! program that applies it afterwards would evaluate one convention's
      ! field against the other's cooling term.
      type(run_options_t), intent(in) :: o

      if (o%ax_field) use_mathis_corrected = .not. o%mathis_orig
      if (o%ax_solver) then
         stoch_method = trim(o%solver)
         qm_method    = trim(o%qm_cooling)
      end if
      if (o%ax_emission) then
         use_induced_emission = o%induced
         gd_photon_cutoff     = o%photon_cutoff
      end if
      if (o%set_nstate) qm_nstate_default = o%nstate
      if (o%set_nisrf)  qm_nisrf_max      = o%nisrf
      if (o%ax_graphite .and. len_trim(o%graphite) > 0) &
         qpah_graphite_source = trim(o%graphite)
      if (o%ax_stage1) s1_density_corrected = o%stage1_density_corrected
      if (o%ax_pah_xsec) qpah_xsec_vintage = trim(o%pah_xsec)
   end subroutine apply_run_options


   subroutine add_hard_euv_component(lam, J_lam)
      ! Replace the band below the Lyman limit -- where the Mathis field is
      ! identically zero -- by a diluted 1e5 K blackbody, so that the field
      ! really does carry photons above 13.6 eV.
      real(wp), intent(in)    :: lam(:)      ! [um]
      real(wp), intent(inout) :: J_lam(:)
      integer :: k
      do k = 1, size(lam)
         if (lam(k) < LAM_LYMAN_UM) J_lam(k) = W_HARD_EUV * bbody(T_HARD_EUV, lam(k))
      end do
   end subroutine add_hard_euv_component

end module sed_apply_options
