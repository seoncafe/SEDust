# SEDust -- Monte Carlo P(T) and SED for stochastically heated grains,
# following Draine & Anderson (1985).
#
# Executables:
#   main_mc.x           single-grain standalone (PIIM Fig. 24.5 / 24.6)
#   main_mc_sed.x       size-distribution SED builder, OpenMP-parallel
#   main_pT_compare.x   single-grain P(T): MC vs the GD matrix solver
#
# All shared modules come from ../sed/src -- this subproject has no
# dependency outside the SEDust tree.

FC       = gfortran
# Compiler family, as in ../tmatrix/Makefile: only the flag spellings differ.
# On the oneAPI releases installed here `ifort` is the IFX compiler.  Set
# FC_FAMILY on the command line for a wrapper whose name hides the compiler.
ifneq (,$(filter ifort ifx mpiifort mpiifx,$(notdir $(FC))))
   FC_FAMILY ?= intel
else
   FC_FAMILY ?= gnu
endif
# Warning set: -Wall -Wextra minus a few false-positive categories.
#   compare-reals: intended sentinel / exactly-0 / NaN comparisons (documented).
#   unused-dummy-argument: kept to preserve shared interfaces.
#   unused-parameter: named constants kept for documentation.
#   do-subscript: false positives from the guarded trapezoid loops.
WARN     = -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument \
           -Wno-unused-parameter -Wno-do-subscript

# ---- HDF5 ---------------------------------------------------------------
# Same switch as ../sed/Makefile: the dust-model modules read the HDF5 optics
# products by default.  make HDF5=0 compiles the HDF5 paths out and the
# readers fall back to the text products.
HDF5        ?= 1
ifeq ($(FC_FAMILY),intel)
   HDF5_PREFIX ?= /data/opt/hdf5_intel
else
   HDF5_PREFIX ?= /data/opt/hdf5_gfortran
endif
ifneq ($(HDF5),0)
   H5FLAGS = -DSEDUST_HDF5 -I$(HDF5_PREFIX)/include
   H5LIBS  = $(HDF5_PREFIX)/lib/libhdf5_fortran.a $(HDF5_PREFIX)/lib/libhdf5.a \
             -lsz -lz -ldl -lm
else
   H5FLAGS =
   H5LIBS  =
endif

FFLAGS   = -O2 -ffree-line-length-none -cpp $(WARN) $(H5FLAGS)
OMPFLAGS = -fopenmp
# Intel spellings (see ../sed/build_lib.sh): -fpp, -qopenmp, -recursive, and
# -fp-model precise in place of Intel's default -fp-model fast.
ifeq ($(FC_FAMILY),intel)
   FFLAGS   = -O2 -fpp -fp-model precise -recursive $(H5FLAGS)
   OMPFLAGS = -qopenmp
endif

SED_DIR = ../sed/src

SHARED_BASE = $(SED_DIR)/constants.f90 $(SED_DIR)/sed_paths.f90 \
              $(SED_DIR)/sed_mathlib.f90
SHARED_GRA  = $(SED_DIR)/q_graphite_d16_sphere.f90
SHARED_RAD  = $(SED_DIR)/radfield.f90
SHARED_ENT  = $(SED_DIR)/enthalpy_v2.f90 $(SED_DIR)/enthalpy_astrodust.f90

# The dust-model modules of ../sed beyond the SHARED_* sets above.
# Keep their relative order the same as SRC_DUST_MODEL in ../sed/Makefile.
SED_DUST_MODEL_SRC = $(SED_DIR)/sedust_h5.f90 $(SED_DIR)/sedust_product.f90 \
                     $(SED_DIR)/size_dist.f90 $(SED_DIR)/q_table.f90 \
                     $(SED_DIR)/q_table_jori.f90 $(SED_DIR)/mie.f90 \
                     $(SED_DIR)/q_astrodust.f90 $(SED_DIR)/q_graphite.f90 \
                     $(SED_DIR)/q_graphite_d16.f90 $(SED_DIR)/qpah.f90 \
                     $(SED_DIR)/p_sub.f90 $(SED_DIR)/stoch_qm.f90 \
                     $(SED_DIR)/pah_ioniz.f90 $(SED_DIR)/grain_dist.f90 \
                     $(SED_DIR)/q_silicate.f90 $(SED_DIR)/pah_ld01.f90 \
                     $(SED_DIR)/scatmat_aligned.f90 $(SED_DIR)/dust_model_mod.f90 \
                     $(SED_DIR)/zubko_io.f90 $(SED_DIR)/kext_table.f90 \
                     $(SED_DIR)/q_component.f90 $(SED_DIR)/dustem_io.f90 \
                     $(SED_DIR)/sed_astrodust.f90

# MC engine + RNG + grain type (shared by all executables)
SRC_ENGINE = src/mc_heatcap.f90 \
             src/mc_rng.f90 \
             src/mc_grain_type.f90 \
             src/mc_engine.f90

SRC_MC          = $(SRC_ENGINE) src/main_mc.f90
SRC_MC_SED      = $(SRC_ENGINE) src/mc_sed.f90 src/main_mc_sed.f90
SRC_PT_COMPARE  = $(SRC_ENGINE) src/main_pT_compare.f90

SRC_MC_ALL         = $(SHARED_BASE) $(SHARED_GRA) $(SHARED_RAD) $(SHARED_ENT) \
                     $(SRC_MC)
SRC_MC_SED_ALL     = $(SHARED_BASE) $(SHARED_GRA) $(SHARED_RAD) $(SHARED_ENT) \
                     $(SED_DUST_MODEL_SRC) $(SRC_MC_SED)
SRC_PT_COMPARE_ALL = $(SHARED_BASE) $(SHARED_GRA) $(SHARED_RAD) $(SHARED_ENT) \
                     $(SED_DUST_MODEL_SRC) $(SRC_PT_COMPARE)

all: main_mc.x main_mc_sed.x main_pT_compare.x

# The drivers write here; git carries no empty directory, so a fresh clone
# would otherwise reach the end of a run and fail on the first open().
output:
	@mkdir -p output

main_mc.x: $(SRC_MC_ALL) | output
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./main_mc.x'

main_mc_sed.x: $(SRC_MC_SED_ALL) | output
	$(FC) $(FFLAGS) $(OMPFLAGS) $^ $(H5LIBS) -o $@
	@echo 'Built: ./main_mc_sed.x (OpenMP)'

main_pT_compare.x: $(SRC_PT_COMPARE_ALL) | output
	$(FC) $(FFLAGS) $(OMPFLAGS) $^ $(H5LIBS) -o $@
	@echo 'Built: ./main_pT_compare.x'

clean:
	rm -f *.o *.mod

cleanall: clean
	rm -f main_mc.x main_mc_sed.x main_pT_compare.x
	rm -f output/Tt_*.dat output/PT_*.dat output/sed_*.dat \
	      output/pT_compare_*.dat

.PHONY: all clean cleanall

# Every executable recipe compiles the shared module sources into this
# directory, so concurrent targets race on the same .mod files (make -j fails
# with .mod0 rename errors). Serialize until the build uses one object/module
# graph in a separate directory.
.NOTPARALLEL:
