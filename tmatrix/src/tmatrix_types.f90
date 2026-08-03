module tmatrix_types
   use, intrinsic :: iso_fortran_env, only: real32, int64
   use tmatrix_kinds, only: wp
   implicit none
   private

   !! Length of every generalized-spherical-function expansion returned by
   !! this library.  It repeats NPL = 2*NPN1+1 of the numerical core, which
   !! is the storage the core fills; the physically populated range is
   !! 1..lmax+1 with lmax = 2*nmax_tm <= 2*NPN4 = 160.
   integer, parameter, public :: tmatrix_expansion_size = 201

   type, public :: tmatrix_options_t
      !! Mishchenko shape convention: shape=-1 is a spheroid and an aspect
      !! ratio greater than one denotes an oblate particle.
      real(wp) :: aspect_ratio = 1.4_wp
      real(wp) :: tolerance    = 1.0e-3_wp
      integer  :: shape        = -1
      integer  :: ndgs         = 2
   end type tmatrix_options_t

   type, public :: tmatrix_scatmat_t
      !! Generalized-spherical-function expansion of the random-orientation
      !! scattering matrix, in the Mishchenko / Hovenier convention with
      !! u = cos(scattering angle) and element (1) the l = 0 term:
      !!
      !!   F11     = sum_l al1(l+1) P_l(u)
      !!   F22+F33 = sum_l (al2+al3)(l+1) P^l_{2,2}(u)
      !!   F22-F33 = sum_l (al2-al3)(l+1) P^l_{2,-2}(u)
      !!   F44     = sum_l al4(l+1) P_l(u)
      !!   F12     = sum_l be1(l+1) P^l_{0,2}(u)
      !!   F34     = sum_l be2(l+1) P^l_{0,2}(u)
      !!
      !! The normalization is al1(1) = 1, i.e. (1/2) integral F11 du = 1
      !! over u in [-1,1].  Entries beyond lmax+1 are zero.
      real(wp) :: al1(tmatrix_expansion_size) = 0.0_wp
      real(wp) :: al2(tmatrix_expansion_size) = 0.0_wp
      real(wp) :: al3(tmatrix_expansion_size) = 0.0_wp
      real(wp) :: al4(tmatrix_expansion_size) = 0.0_wp
      real(wp) :: be1(tmatrix_expansion_size) = 0.0_wp
      real(wp) :: be2(tmatrix_expansion_size) = 0.0_wp
      !! Truncation order of the expansion; the highest populated index is
      !! lmax+1.
      integer  :: lmax    = 0
      !! Multipole truncation order of the converged T-matrix left in the
      !! workspace.  It is the highest N, and the highest azimuthal order M,
      !! for which full_tstore is filled, hence the NMAX a fixed-orientation
      !! amplitude evaluation must use.
      integer  :: nmax_tm = 0
      !! Generation stamp of the T-matrix this expansion was obtained with: a
      !! copy of work%tm_tag as the evaluation left it.  A fixed-orientation
      !! amplitude evaluation compares it against the workspace, so a T-matrix
      !! replaced by a later evaluation is reported instead of read.  The
      !! default 0 never matches a workspace that has evaluated anything.
      integer(int64) :: state_tag = 0
   end type tmatrix_scatmat_t

   type, public :: tmatrix_result_t
      real(wp) :: qext        = 0.0_wp
      real(wp) :: qsca        = 0.0_wp
      real(wp) :: qabs        = 0.0_wp
      real(wp) :: albedo      = 0.0_wp
      real(wp) :: asymmetry   = 0.0_wp
      integer  :: status      = 0
      integer  :: legacy_ierr = 0
      integer  :: lapack_info = 0
      character(len=160) :: message = ''
   end type tmatrix_result_t

   type, public :: tmatrix_workspace_t
      !! Caller-owned mutable state for the default reentrant backend.  CT,
      !! CTT, Bessel, T-matrix, and GSP storage are distinct components; no
      !! legacy overlay is reconstructed.  Different workspaces may be used
      !! concurrently, but a single workspace must have one active caller.
      logical :: initialized = .false.
      ! Identity of the T-matrix currently held in full_tstore.  tm_tag counts
      ! evaluations on this workspace and advances on a failed one too, because
      ! a failed evaluation may already have overwritten part of the storage;
      ! tm_filled says whether the last evaluation left a converged T-matrix;
      ! tm_lambda_um and tm_nmax are the wavelength and the multipole
      ! truncation order it was computed at.  A fixed-orientation amplitude
      ! evaluation reads all four instead of trusting caller-supplied values.
      integer(int64) :: tm_tag       = 0
      logical        :: tm_filled    = .false.
      real(wp)       :: tm_lambda_um = 0.0_wp
      integer        :: tm_nmax      = 0
      real(wp) :: fac(900) = 0.0_wp
      real(wp) :: ssign(900) = 0.0_wp
      real(wp), allocatable :: ct_tr(:,:), ct_ti(:,:)
      real(wp), allocatable :: ctt_qr(:,:), ctt_qi(:,:)
      real(wp), allocatable :: ctt_rgqr(:,:), ctt_rgqi(:,:)
      real(wp), allocatable :: cbess_j(:,:), cbess_y(:,:), cbess_jr(:,:), cbess_ji(:,:)
      real(wp), allocatable :: cbess_dj(:,:), cbess_dy(:,:), cbess_djr(:,:), cbess_dji(:,:)
      real(real32), allocatable :: cbess_gsp_b1r(:,:,:), cbess_gsp_b1i(:,:,:)
      real(real32), allocatable :: cbess_gsp_b2r(:,:,:), cbess_gsp_b2i(:,:,:)
      ! Default full-direct state.  Deliberately distinct T and GSP arrays
      ! replace the historical incompatible overlay.  Because the expansion
      ! routine only reads full_tstore, the converged T-matrix survives an
      ! evaluation intact and a fixed-orientation amplitude evaluation can
      ! read it afterwards.  Slot k of the fourth dimension holds, in order,
      ! Re T11, Re T12, Re T21, Re T22, Im T11, Im T12, Im T21, Im T22.
      real(real32), allocatable :: full_tstore(:,:,:,:)
      real(wp), allocatable :: full_tmatr_work(:,:,:)
      real(real32), allocatable :: full_gsp_d1(:,:,:), full_gsp_d2(:,:,:)
      real(real32), allocatable :: full_gsp_d3(:,:,:), full_gsp_d4(:,:,:)
      real(real32), allocatable :: full_gsp_d5r(:,:,:), full_gsp_d5i(:,:,:)
      ! Largest scratch arrays of TT, TMATR0/TMATR, and GSP.  They were
      ! automatic arrays living on the caller's stack, which forced an OpenMP
      ! worker stack of several megabytes; owning them here leaves the deepest
      ! full-direct call chain well under half a megabyte of stack.  TMATR0 and
      ! TMATR are never active at the same time and share one D1/D2 pair.
      complex(wp), allocatable :: tt_zq(:,:)
      real(wp), allocatable :: tmatr_d1(:,:), tmatr_d2(:,:)
      real(wp), allocatable :: gsp_tr1(:,:), gsp_tr2(:,:), gsp_ti1(:,:), gsp_ti2(:,:)
      real(wp), allocatable :: gsp_g1(:,:), gsp_g2(:,:)
      real(wp), allocatable :: gsp_fr(:,:), gsp_fi(:,:), gsp_ff(:,:)
   end type tmatrix_workspace_t
end module tmatrix_types
