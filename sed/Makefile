# SEDust -- SED solver build.
#
# Targets:
#   make                     main_astrodust.x main_dl07.x main_zubko.x
#                            make_enthalpy.x calc_kext.x
#   make libsedust.a         the model-agnostic dust-emission library
#   make calc_kext.x         size-integrated transport optics for any model
#   make clean / cleanall
#
# Run (from this directory; data paths are relative to ../):
#   ./make_enthalpy.x
#   ./main_astrodust.x                 # GD solver (default)
#   ./main_astrodust.x qm              # energy-space, thermal-discrete
#   ./main_astrodust.x qm_dbcon        # energy-space, thermal-continuous
#   ./main_dl07.x
#   ./main_zubko.x                     # Zubko/ZDA emission
#   ./main_zubko.x euv                 # ... with the field carried into the EUV
#   ./calc_kext.x astrodust            # -> ../data/kext_astrodust_MW.dat
#   ./calc_kext.x astrodust euv        # carried into the ionizing band
#   ./calc_kext.x dl07 euv
#   ./calc_kext.x zubko
#   ./calc_kext.x from_files <descriptor> [data_dir]

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

# ---- T-matrix library ------------------------------------------------
# The astrodust model's extreme-ultraviolet band (a lam_min shorter than the Q
# table's 0.0912 um) is solved with the T-matrix, on the same b/a = 1.4 oblate
# spheroid the table itself is, so everything that links sed_astrodust.f90 now
# also links libtmatrix.a.  Build it first:
#     cd ../tmatrix && make libtmatrix.a
TMDIR   = ../tmatrix
TMARCH  = $(TMDIR)/libtmatrix.a
TMLIB   = -L$(TMDIR) -ltmatrix
FFLAGS += -I$(TMDIR)

# The two driver modules that own the size-parameter regime selection
# (Rayleigh dipole / T-matrix / geometric optics).  They live outside the
# archive, are compiled here, and must precede src/sed_astrodust.f90.
SRC_TMDRV = $(TMDIR)/driver/asymptotic_optics.f90 $(TMDIR)/driver/spheroid_optics.f90

# ---- source groups, in dependency order ------------------------------
SRC_BASE    = src/constants.f90 src/sed_mathlib.f90

# DL07-support modules (also used by sed_astrodust.f90, so required by every
# target that links it): PAH ionization (WD01b charging), WD01/DL07 size
# distribution, silicate optics.  Must precede sed_astrodust.f90.
SRC_DL07MOD = src/pah_ioniz.f90 src/grain_dist.f90 src/q_silicate.f90 \
              src/pah_ld01.f90

# Everything the SED pipeline needs, minus the driver.
SRC_PIPE    = $(SRC_BASE) src/enthalpy_v2.f90 \
              src/size_dist.f90 src/q_table.f90 \
              src/enthalpy_astrodust.f90 \
              src/mie.f90 src/q_astrodust.f90 \
              src/q_graphite.f90 src/q_graphite_d16.f90 \
              src/q_graphite_d16_sphere.f90 src/qpah.f90 \
              src/radfield.f90 src/p_sub.f90 \
              src/stoch_qm.f90 $(SRC_DL07MOD) \
              src/dust_model_mod.f90 src/kext_table.f90 src/zubko_io.f90 \
              $(SRC_TMDRV) src/sed_astrodust.f90

SRC_ENTH    = $(SRC_BASE) src/enthalpy_v2.f90 \
              src/enthalpy_astrodust.f90 src/make_enthalpy.f90

all: make_enthalpy.x main_astrodust.x main_dl07.x main_zubko.x calc_kext.x

# Fail with an instruction rather than a link error when the T-matrix archive is
# absent.  No prerequisites, so this recipe runs only if the file does not
# exist.  It sits after `all` so that `all` stays the default goal.
$(TMARCH):
	@echo '*** $(TMARCH) not found.'
	@echo '*** Build it first:   cd $(TMDIR) && make libtmatrix.a'
	@false

make_enthalpy.x: $(SRC_ENTH)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./make_enthalpy.x'

main_astrodust.x: $(SRC_PIPE) src/main_astrodust.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $^ $(TMLIB) -o $@
	@echo 'Built: ./main_astrodust.x'

main_dl07.x: $(SRC_PIPE) src/main_dl07.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $^ $(TMLIB) -o $@
	@echo 'Built: ./main_dl07.x'

# Zubko/ZDA emission. The only shipped model whose optics grid reaches far
# past the illuminated band, so it is the one that exercises the difference
# between "hardest photon the field carries" and "short end of the grid".
main_zubko.x: $(SRC_PIPE) src/main_zubko.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $^ $(TMLIB) -o $@
	@echo 'Built: ./main_zubko.x'

# Size-integrated transport optics (lambda, albedo, <cos>, C_ext/C_abs/C_sca
# per H) for astrodust, DL07, Zubko or any file-defined model -> ../data/.
# Reuses the full SED pipeline, so every table it writes is the same optics
# the RT hosts get from dust_extinction; only the driver differs.
calc_kext.x: $(SRC_PIPE) src/calc_kext.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $^ $(TMLIB) -o $@
	@echo 'Built: ./calc_kext.x'

# ---- libsedust.a: the model-agnostic dust-emission library, for a Fortran
# 3D RT code to link.  Usage from the RT:  use dust_lib ;  link
#   -L<sed> -lsedust -L<tmatrix> -ltmatrix
# with the .mod search paths  -I<sed> -I<tmatrix>.  BOTH archives are needed:
# the astrodust EUV band calls the T-matrix.
LIB_SRC = $(SRC_PIPE) src/dust_lib.f90
LIB_OBJ = $(notdir $(LIB_SRC:.f90=.o))

libsedust.a: $(LIB_SRC) | $(TMARCH)
	$(FC) $(FFLAGS) -c $(LIB_SRC)
	rm -f $@
	ar rcs $@ $(LIB_OBJ)
	@echo 'Built: libsedust.a  (link with -L. -lsedust -L$(TMDIR) -ltmatrix -I. -I$(TMDIR))'

clean:
	rm -f *.o *.mod

cleanall: clean
	rm -f *.a make_enthalpy.x main_astrodust.x main_dl07.x main_zubko.x calc_kext.x
	rm -f output/astrodust_irem_ours_*.dat output/enthalpy_*.dat

.PHONY: all clean cleanall

# Every executable recipe compiles the shared module sources into this
# directory, so concurrent targets race on the same .mod files (make -j fails
# with .mod0 rename errors). Serialize until the build uses one object/module
# graph in a separate directory.
.NOTPARALLEL:
