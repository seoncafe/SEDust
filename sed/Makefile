# SEDust -- SED solver build.
#
# Targets:
#   make                     main_astrodust.x main_dl07.x make_enthalpy.x
#   make libsedust.a         the model-agnostic dust-emission library
#   make calc_kext_astrodust.x   astrodust+PAH extinction/albedo/<cos>/K_abs
#   make calc_kext_dl07.x        DL07 (WD01 MW) extinction/albedo/<cos>
#   make calc_polext.x           HD23 polarized extinction per H, vs release
#   make clean / cleanall
#
# Run (from this directory; data paths are relative to ../):
#   ./make_enthalpy.x
#   ./main_astrodust.x                 # GD solver (default)
#   ./main_astrodust.x qm              # energy-space, thermal-discrete
#   ./main_astrodust.x qm_dbcon        # energy-space, thermal-continuous
#   ./main_dl07.x

FC      = gfortran
# Portable baseline by default so libsedust.a runs on any node it is copied
# to.  Uncomment the -march=native line for a host-optimized build (results
# can differ in the last digits through FMA contraction).
# Warning set: -Wall -Wextra minus a few false-positive categories.
#   compare-reals: intended sentinel / exactly-0 / NaN comparisons (documented).
#   unused-dummy-argument: kept to preserve shared interfaces.
#   unused-parameter: named constants kept for documentation.
#   do-subscript: false positives from the guarded trapezoid loops.
WARN    = -Wall -Wextra -Wno-compare-reals -Wno-unused-dummy-argument \
          -Wno-unused-parameter -Wno-do-subscript
FFLAGS  = -O3 -ffree-line-length-none $(WARN)
# FFLAGS  = -O3 -march=native -ffree-line-length-none $(WARN)
FFLAGS += -fopenmp

# ---- source groups, in dependency order ------------------------------
SRC_BASE    = src/constants.f90 src/sed_mathlib.f90

# DL07-support modules (also used by sed_astrodust.f90, so required by every
# target that links it): PAH ionization (WD01b charging), WD01/DL07 size
# distribution, silicate optics.  Must precede sed_astrodust.f90.
SRC_DL07MOD = src/pah_ioniz.f90 src/grain_dist.f90 src/q_silicate.f90 \
              src/pah_ld01.f90

# Everything the SED pipeline needs, minus the driver.
SRC_PIPE    = $(SRC_BASE) src/enthalpy_v2.f90 \
              src/size_dist.f90 src/q_table.f90 src/q_table_jori.f90 \
              src/enthalpy_astrodust.f90 \
              src/mie.f90 src/q_graphite.f90 src/q_graphite_d16.f90 \
              src/q_graphite_d16_sphere.f90 src/qpah.f90 \
              src/radfield.f90 src/p_sub.f90 \
              src/stoch_qm.f90 $(SRC_DL07MOD) \
              src/dust_model_mod.f90 src/zubko_io.f90 \
              src/sed_astrodust.f90

SRC_ENTH    = $(SRC_BASE) src/enthalpy_v2.f90 \
              src/enthalpy_astrodust.f90 src/make_enthalpy.f90

# Cross section / albedo / <cos> for the DL07 (WD01 MW) model.  Reuses the
# same optics the SED model uses (silicate + graphite Mie, qpah xi-blend).
SRC_KEXT    = $(SRC_BASE) src/mie.f90 \
              src/q_graphite.f90 src/q_graphite_d16.f90 \
              src/q_graphite_d16_sphere.f90 \
              src/qpah.f90 src/pah_ioniz.f90 src/grain_dist.f90 \
              src/q_silicate.f90 src/calc_kext_dl07.f90

all: make_enthalpy.x main_astrodust.x main_dl07.x

make_enthalpy.x: $(SRC_ENTH)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./make_enthalpy.x'

main_astrodust.x: $(SRC_PIPE) src/main_astrodust.f90
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./main_astrodust.x'

main_dl07.x: $(SRC_PIPE) src/main_dl07.f90
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./main_dl07.x'

calc_kext_dl07.x: $(SRC_KEXT)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./calc_kext_dl07.x'

# astrodust+PAH extinction/albedo/<cos>/K_abs -> ../data/kext_astrodust_MW.dat.
# Reuses the full SED pipeline (sed_init) so the optics are identical to
# main_astrodust.x; just swaps the driver.
calc_kext_astrodust.x: $(SRC_PIPE) src/calc_kext_astrodust.f90
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./calc_kext_astrodust.x'

# Polarized extinction per H from the orientation-resolved DH21 spheroid
# table, compared against ../data/release/polarized_extinction.dat.
SRC_POLEXT = $(SRC_BASE) src/size_dist.f90 src/q_table_jori.f90 \
             src/calc_polext.f90

calc_polext.x: $(SRC_POLEXT)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./calc_polext.x'

# ---- libsedust.a: the model-agnostic dust-emission library, for a Fortran
# 3D RT code to link.  Usage from the RT:  use dust_lib ;  link
#   -L<sed> -lsedust  with the .mod search path  -I<sed>
LIB_SRC = $(SRC_PIPE) src/dust_lib.f90
LIB_OBJ = $(notdir $(LIB_SRC:.f90=.o))

libsedust.a: $(LIB_SRC)
	$(FC) $(FFLAGS) -c $(LIB_SRC)
	rm -f $@
	ar rcs $@ $(LIB_OBJ)
	@echo 'Built: libsedust.a  (link with -L. -lsedust -I.)'

clean:
	rm -f *.o *.mod

cleanall: clean
	rm -f *.a make_enthalpy.x main_astrodust.x main_dl07.x \
	      calc_kext_dl07.x calc_kext_astrodust.x calc_polext.x
	rm -f output/astrodust_irem_ours_*.dat output/enthalpy_*.dat

.PHONY: all clean cleanall
