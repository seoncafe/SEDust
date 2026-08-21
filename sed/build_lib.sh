#!/bin/bash
# Build libsedust.a, the model-agnostic dust-emission library, for a Fortran
# 3D radiative-transfer host to link.  This is the script an RT tree copies
# next to its own copy of sed/src/; the Makefile in this directory builds the
# same library plus the CLI drivers.
#
# Two builds, because only ONE thing in the library needs the T-matrix -- the
# astrodust extreme-ultraviolet band, i.e. a lam_min shorter than the Q table's
# 0.0912 um (13.6 eV) end, solved on the same oblate spheroid the table is:
#
#   ./build_lib.sh                 No T-matrix anywhere in the archive.  A host
#                                  that transports no ionizing radiation, or
#                                  that is content with the volume-equivalent
#                                  sphere (euv_tmatrix = .false.), needs
#                                  nothing else.  euv_tmatrix = .true. is then
#                                  REFUSED -- status 6 from build_astrodust, 8
#                                  from build_dust, which keeps its own single
#                                  vocabulary -- never answered with the sphere
#                                  behind the caller's back.
#
#   WITH_TMATRIX=1 ./build_lib.sh  The T-matrix and the spheroid optics of the
#                                  EUV band go into the SAME archive, so the
#                                  host still links only -lsedust.  It must
#                                  then call
#                                      use euv_astrodust_tmatrix, only: &
#                                          use_tmatrix_euv_band_optics
#                                      call use_tmatrix_euv_band_optics()
#                                  once before build_astrodust.
#
# Output: sed/lib/libsedust.a and the .mod files the host compiles against
# (one search path, -Ilib).  Re-run after updating sed/src/.
#
# Overridable from the environment: FC, FLAGS, MODOUT, TMDIR, and the T-matrix
# flag sets TM_FREE / TM_CORE / TM_FIXED.  The defaults below are ifort's; a
# gfortran FC selects gfortran's.
set -e
cd "$(dirname "$0")"                 # sed/

FC=${FC:-ifort}
TMDIR=${TMDIR:-../tmatrix}

# ---- HDF5 ------------------------------------------------------------
# build_dust reads the optics products out of data/<model>/sedust_<model>.h5, so the
# archive carries the HDF5 layer.  ON by default, and switched off for a host
# that has no HDF5 or wants a text-only archive:
#
#   HDF5=0 ./build_lib.sh          the HDF5 paths compile out; build_dust then
#                                  falls back to the text products
#
# A host linking the default archive must put -lsedust BEFORE the HDF5
# libraries on its link line: static archives resolve left to right, and it is
# libsedust.a that needs the HDF5 symbols.
HDF5=${HDF5:-1}
case "$(basename "$FC")" in
   gfortran*|gfort*)  HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_gfortran} ;;
   *)                 HDF5_PREFIX=${HDF5_PREFIX:-/data/opt/hdf5_intel} ;;
esac
if [ "$HDF5" != "0" ]; then
   H5FLAGS="-DSEDUST_HDF5 -I$HDF5_PREFIX/include"
else
   H5FLAGS=""
fi

case "$(basename "$FC")" in
   gfortran*|gfort*)
      : "${MODOUT:=-Jlib}"
      # -frecursive: the EUV band loop in sed_astrodust.f90 is threaded and
      # reaches the T-matrix through spheroid_optics, so local variables belong
      # on the stack.  Enabling OpenMP already implies this, on both compilers;
      # stating it keeps a threaded path from resting on that side effect.  It
      # costs nothing here -- GFortran stacks local arrays either way.
      : "${FLAGS:=-O3 -ffree-line-length-none -fopenmp -frecursive -cpp -w}"
      # The T-matrix sources carry their own flags, matching tmatrix/Makefile:
      # -frecursive keeps every migrated routine reentrant (a workspace belongs
      # to one caller, and callers run concurrently), and -fno-inline on
      # tmatrix_core fixes where the compiler contracts multiply-add pairs,
      # which is what makes its convergence test reproducible.  The archive is
      # deliberately built without OpenMP: the threading is on the SEDust side.
      : "${TM_FREE:=-O2 -ffree-line-length-none -frecursive -w}"
      : "${TM_CORE:=$TM_FREE -fno-inline}"
      : "${TM_FIXED:=-O2 -ffixed-line-length-none -fmax-stack-var-size=65536 -frecursive -w}"
      ;;
   *)
      : "${MODOUT:=-module lib}"
      # -recursive for the reason given in the GNU branch above; on Intel it is
      # not free, because the default (-auto-scalar) leaves fixed-size local
      # arrays in static storage.  Measured: an object holding one 4096-element
      # local array carries 32 kB of .bss without it and none with it.
      #
      # -fp-model precise, and no -xHost.  Intel defaults to -fp-model fast,
      # which reassociates and contracts freely.  The T-matrix convergence test
      # compares successive iterates against DDELT, so value-unsafe
      # rearrangement moves the iteration count: at the default the 120-case
      # golden regression misses its 1e-10 tolerance and both golden tests fail.
      # The SED sources do not run that test themselves, but they are compiled
      # into one archive with the routines that do, and euv_astrodust_tmatrix
      # calls straight into them.
      : "${FLAGS:=-O3 -qopenmp -recursive -fp-model precise -fpp -w}"
      : "${TM_FREE:=-O2 -recursive -fp-model precise -w}"
      # Intel counterpart of -fno-inline.  The reproducibility argument recorded
      # in the GNU branch is a GNU-to-GNU one and does not carry across
      # compilers; the flag keeps the Intel build's core un-inlined all the same.
      : "${TM_CORE:=$TM_FREE -inline-level=0}"
      # lpd.f is fixed-form and runs past column 72.
      : "${TM_FIXED:=-O2 -extend-source 132 -recursive -fp-model precise -w}"
      ;;
esac

# SEDust sources, in dependency order.  `sed_mathlib` is named so upstream to
# avoid a clash with an RT host's own `mathlib` module.
SED_SRCS="constants sed_paths sed_mathlib sedust_h5 sedust_product enthalpy_v2 size_dist q_table \
      enthalpy_astrodust mie \
      q_astrodust q_graphite q_graphite_d16 q_graphite_d16_sphere qpah radfield \
      p_sub stoch_qm pah_ioniz grain_dist q_silicate pah_ld01 dust_model_mod \
      kext_table zubko_io q_component sed_astrodust dust_lib"

mkdir -p lib
rm -f lib/*.o lib/*.mod lib/libsedust.a

if [ -n "$WITH_TMATRIX" ]; then
   # T-matrix first: euv_astrodust_tmatrix at the end needs both these modules
   # and sed_astrodust_mod.
   for f in tmatrix_kinds tmatrix_status tmatrix_types; do
      echo "  FC $f"
      $FC $TM_FREE $MODOUT -Ilib -c $TMDIR/src/$f.f90 -o lib/$f.o
   done
   echo "  FC tmatrix_core"
   $FC $TM_CORE $MODOUT -Ilib -c $TMDIR/src/tmatrix_core.f90 -o lib/tmatrix_core.o
   for f in tmatrix_full_direct_bridge fixed_orientation_amplitude \
            fixed_orientation_phase_matrix scattering_matrix_expansion tmatrix_api; do
      echo "  FC $f"
      $FC $TM_FREE $MODOUT -Ilib -c $TMDIR/src/$f.f90 -o lib/$f.o
   done
   echo "  FC lpd"
   $FC $TM_FIXED -c $TMDIR/src/lpd.f -o lib/lpd.o
   # The two modules that own the size-parameter regime selection (Rayleigh
   # dipole / T-matrix / geometric optics).
   for f in asymptotic_optics spheroid_optics; do
      echo "  FC $f"
      $FC $TM_FREE $MODOUT -Ilib -c $TMDIR/driver/$f.f90 -o lib/$f.o
   done
fi

for f in $SED_SRCS; do
   echo "  FC $f"
   $FC $FLAGS $H5FLAGS $MODOUT -Ilib -c src/$f.f90 -o lib/$f.o
done

if [ -n "$WITH_TMATRIX" ]; then
   echo "  FC euv_astrodust_tmatrix"
   $FC $FLAGS $MODOUT -Ilib -c src/euv_astrodust_tmatrix.f90 -o lib/euv_astrodust_tmatrix.o
fi

ar rcs lib/libsedust.a lib/*.o
if [ -n "$WITH_TMATRIX" ]; then
   echo "Built: sed/lib/libsedust.a  (with the T-matrix; call use_tmatrix_euv_band_optics)"
else
   echo "Built: sed/lib/libsedust.a  (no T-matrix; euv_tmatrix = .true. gives"
   echo "       status 6 from build_astrodust, 8 from build_dust)"
fi
