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
# Warning set: -Wall -Wextra minus a few false-positive categories.
#   compare-reals: intended sentinel / exactly-0 / NaN comparisons (documented).
#   unused-dummy-argument: kept to preserve shared interfaces.
#   unused-parameter: named constants kept for documentation.
#   do-subscript: false positives from the guarded trapezoid loops.
WARN     = -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument \
           -Wno-unused-parameter -Wno-do-subscript
FFLAGS   = -O2 -ffree-line-length-none $(WARN)
OMPFLAGS = -fopenmp

SED_DIR = ../sed/src

SHARED_BASE = $(SED_DIR)/constants.f90 $(SED_DIR)/sed_mathlib.f90
SHARED_GRA  = $(SED_DIR)/q_graphite_d16_sphere.f90
SHARED_RAD  = $(SED_DIR)/radfield.f90
SHARED_ENT  = $(SED_DIR)/enthalpy_v2.f90 $(SED_DIR)/enthalpy_astrodust.f90

# Astrodust pipeline modules pulled in for the SED comparison
SED_PIPELINE = $(SED_DIR)/size_dist.f90 $(SED_DIR)/q_table.f90 \
               $(SED_DIR)/q_table_jori.f90 \
               $(SED_DIR)/mie.f90 $(SED_DIR)/q_graphite.f90 \
               $(SED_DIR)/q_graphite_d16.f90 $(SED_DIR)/qpah.f90 \
               $(SED_DIR)/p_sub.f90 $(SED_DIR)/stoch_qm.f90 \
               $(SED_DIR)/pah_ioniz.f90 $(SED_DIR)/grain_dist.f90 \
               $(SED_DIR)/q_silicate.f90 $(SED_DIR)/pah_ld01.f90 \
               $(SED_DIR)/scatmat_aligned.f90 \
               $(SED_DIR)/dust_model_mod.f90 $(SED_DIR)/zubko_io.f90 \
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
                     $(SED_PIPELINE) $(SRC_MC_SED)
SRC_PT_COMPARE_ALL = $(SHARED_BASE) $(SHARED_GRA) $(SHARED_RAD) $(SHARED_ENT) \
                     $(SED_PIPELINE) $(SRC_PT_COMPARE)

all: main_mc.x main_mc_sed.x main_pT_compare.x

# The drivers write here; git carries no empty directory, so a fresh clone
# would otherwise reach the end of a run and fail on the first open().
output:
	@mkdir -p output

main_mc.x: $(SRC_MC_ALL) | output
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./main_mc.x'

main_mc_sed.x: $(SRC_MC_SED_ALL) | output
	$(FC) $(FFLAGS) $(OMPFLAGS) $^ -o $@
	@echo 'Built: ./main_mc_sed.x (OpenMP)'

main_pT_compare.x: $(SRC_PT_COMPARE_ALL) | output
	$(FC) $(FFLAGS) $(OMPFLAGS) $^ -o $@
	@echo 'Built: ./main_pT_compare.x'

clean:
	rm -f *.o *.mod

cleanall: clean
	rm -f main_mc.x main_mc_sed.x main_pT_compare.x
	rm -f output/Tt_*.dat output/PT_*.dat output/sed_*.dat \
	      output/mc_sed_*.dat output/pT_compare_*.dat

.PHONY: all clean cleanall
