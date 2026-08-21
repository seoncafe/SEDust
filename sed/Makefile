# SEDust -- SED solver build.
#
# Targets:
#   make                     calc_sed.x calc_enthalpy.x calc_kext.x
#   make libsedust.a         the model-agnostic dust-emission library, with no
#                            T-matrix in it (see below)
#   WITH_TMATRIX=1 make libsedust.a
#                            ... with the T-matrix and the spheroid optics of
#                            the astrodust EUV band in the same archive
#   make calc_kext.x         size-integrated transport optics for any model
#   make clean / cleanall
#
# Run (from this directory; data paths are relative to ../):
#   ./calc_enthalpy.x            # astrodust enthalpy tables (add c2 for the
#                                #   Stage-1 density-corrected prefactor)
#   ./calc_qtable.x                    # -> ../data/<model>/ (every model's Q tables)
#   ./calc_sed.x astrodust             # astrodust+PAH SED  -> output/
#   ./calc_sed.x dl07                  # Draine & Li (2007) SED
#   ./calc_sed.x mrn                   # Mathis, Rumpl & Nordsieck (1977) SED
#   ./calc_sed.x zubko                 # Zubko/ZDA SED
#
#   Every program here takes ONE vocabulary of run settings
#   (src/sed_run_options.f90) and declares which axes it has a referent for;
#   each setting tags the output filename.  calc_sed.x declares:
#     solver    heuristic (default) | draine | equil | qm | qm_dbcon | qm_stati
#     grid      euv                  solve the ionizing band as well
#     field     mathis_orig, logU=X, hardfield
#     emission  induced, photcut
#     qm sizes  nstate=N, nisrf=N
#     graphite  gra_d03_sphere | gra_d16_sphere | gra_d16_spheroid
#               (astrodust, dl07: the PAH xi blend)
#     normaliz. dl84 (default) | mrn77   (mrn: which published a^-3.5
#               normalization)
#     enthalpy  c2                   (astrodust: Stage-1 prefactor)
#   ./calc_sed.x astrodust qm_stati nstate=500
#   ./calc_sed.x dl07 lmc2_10 draine euv
#   ./calc_sed.x mrn draine euv
#   ./calc_sed.x zubko euv hardfield
#   ./calc_kext.x astrodust            # -> ../data/kext_astrodust_MW.dat
#   ./calc_kext.x astrodust euv        # carried into the ionizing band
#   ./calc_kext.x dl07 euv
#   ./calc_kext.x mrn
#   ./calc_kext.x mrn mrn77 euv
#   ./calc_kext.x zubko
#   ./calc_kext.x from_files <descriptor> [data_dir]
#   ./check_build_dust.x               # build_dust on HDF5 vs the text builders

FC      = gfortran

# ---- HDF5 -------------------------------------------------------------
# The optics products (kext curves, Q tables, and in v1.20 the polarized
# matrices) are written and read as HDF5; the text products are still written
# beside them, for reading with an eye and for diff.  HDF5 is therefore ON by
# default and can be switched off for a build that has no library:
#
#   make HDF5=0                      text only; the HDF5 paths compile out
#   make HDF5_PREFIX=/somewhere      a different installation
#
# The prefix must hold include/hdf5.mod built by THIS compiler -- a .mod from
# another gfortran release is not readable, which is the usual failure here.
HDF5        ?= 1
ifeq ($(FC),ifort)
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
# The HDF5 layer and the reader that knows the SEDust layout inside a file.
# Both are always compiled: without -DSEDUST_HDF5 every entry point reports
# "not available" and the callers fall through to the text products, so no
# consumer carries an #ifdef of its own.
SRC_H5 = src/sedust_h5.f90 src/sedust_product.f90
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
FFLAGS  = -O3 -ffree-line-length-none -cpp $(WARN) $(H5FLAGS)
# FFLAGS  = -O3 -march=native -ffree-line-length-none $(WARN)
FFLAGS += -fopenmp
# Local variables belong on the run-time stack, not in static storage: the EUV
# band loop in sed_astrodust.f90 is threaded, and it reaches spheroid_optics.f90
# and asymptotic_optics.f90, which this Makefile recompiles with its own flags.
# The flag is a no-op for GFortran, which already stacks local arrays (measured:
# an object holding a 4096-element local array has an empty .bss either way).
# It is not a no-op for Intel, where the same object carries 32 kB of .bss
# without it.  Enabling OpenMP happens to imply it on both compilers, so nothing
# is broken today; stating it keeps the reentrancy of a threaded path from
# resting on that side effect.
FFLAGS += -frecursive

# ---- T-matrix library, for the astrodust EUV band only ---------------
# The astrodust model's extreme-ultraviolet band (a lam_min shorter than the Q
# table's 0.0912 um) is solved on the same oblate spheroid the table itself is,
# which is a T-matrix calculation.  It is the only thing in the tree that needs
# libtmatrix.a, and src/sed_astrodust.f90 reaches it through a procedure
# pointer rather than a use statement, so:
#
#   * `make libsedust.a` links and archives NO T-matrix at all.  A model built
#     from that library answers euv_tmatrix = .true. with status 6.
#   * `WITH_TMATRIX=1 make libsedust.a` adds the T-matrix and the spheroid
#     route to the same archive, so an RT host still links only -lsedust.  The
#     host calls use_tmatrix_euv_band_optics() once before build_astrodust.
#   * the two astrodust CLI programs (calc_sed.x, calc_kext.x) always
#     take the spheroid route, so they are always built against
#     $(TMDIR)/libtmatrix.a.  Build it first:
#         cd ../tmatrix && make libtmatrix.a
TMDIR   = ../tmatrix
TMARCH  = $(TMDIR)/libtmatrix.a
TMINC   = -I$(TMDIR)
TMLIB   = -L$(TMDIR) -ltmatrix

# The spheroid route: the two driver modules that own the size-parameter regime
# selection (Rayleigh dipole / T-matrix / geometric optics), and the module
# that solves the EUV band with them.  The last one uses sed_astrodust_mod, so
# it must FOLLOW src/sed_astrodust.f90.
SRC_EUVTM = $(TMDIR)/driver/asymptotic_optics.f90 \
            $(TMDIR)/driver/spheroid_optics.f90 \
            src/euv_astrodust_tmatrix.f90

# ---- source groups, in dependency order ------------------------------
SRC_BASE    = src/constants.f90 src/sed_paths.f90 src/sed_mathlib.f90

# DL07-support modules (also used by sed_astrodust.f90, so required by every
# target that links it): PAH ionization (WD01b charging), WD01/DL07 size
# distribution, silicate optics.  Must precede sed_astrodust.f90.
SRC_DL07MOD = src/pah_ioniz.f90 src/grain_dist.f90 src/q_silicate.f90 \
              src/pah_ld01.f90

# Everything the SED pipeline needs, minus the driver.  SRC_H5 sits after
# SRC_BASE because src/sedust_product.f90 uses constants, and before
# src/q_table.f90, which reads the astrodust table through it.
SRC_PIPE    = $(SRC_BASE) $(SRC_H5) src/enthalpy_v2.f90 \
              src/size_dist.f90 src/q_table.f90 \
              src/enthalpy_astrodust.f90 \
              src/mie.f90 src/q_astrodust.f90 \
              src/q_graphite.f90 src/q_graphite_d16.f90 \
              src/q_graphite_d16_sphere.f90 src/qpah.f90 \
              src/radfield.f90 src/p_sub.f90 \
              src/stoch_qm.f90 $(SRC_DL07MOD) \
              src/dust_model_mod.f90 src/kext_table.f90 src/zubko_io.f90 \
              src/q_component.f90 \
              src/sed_astrodust.f90

SRC_ENTH    = $(SRC_BASE) src/enthalpy_v2.f90 \
              src/enthalpy_astrodust.f90 $(SRC_CLIOPT) src/calc_enthalpy.f90

# The command-line SYSTEM every program here shares: one vocabulary, one
# parser, one refusal, one filename-tag rule, with each program declaring which
# axes it has a referent for.  It reads the module state of the solver, the
# radiation field and the optics, so it follows SRC_PIPE wherever the pipeline
# is linked; it is NOT part of libsedust.a, a host having no command line.
SRC_CLIOPT  = src/sed_run_options.f90
# Where a setting LANDS: the module state the solver, the field and the optics
# read.  Only the programs that HAVE those modules link it.
SRC_CLIAPPLY = src/sed_apply_options.f90

all: calc_enthalpy.x calc_sed.x calc_kext.x

# Fail with an instruction rather than a link error when the T-matrix archive is
# absent.  No prerequisites, so this recipe runs only if the file does not
# exist.  It sits after `all` so that `all` stays the default goal.
$(TMARCH):
	@echo '*** $(TMARCH) not found.'
	@echo '*** Build it first:   cd $(TMDIR) && make libtmatrix.a'
	@false

calc_enthalpy.x: $(SRC_ENTH)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./calc_enthalpy.x'

# Stored (lambda, a_eff) cross-section tables for the DL07 and Zubko models,
# the counterpart of tmatrix/output/q_astrodust_*.dat.  Both grids of each
# model, computed by this tree's own optics routines.
calc_qtable.x: $(SRC_PIPE) $(SRC_CLIOPT) src/calc_qtable.f90
	$(FC) $(FFLAGS) $^ $(H5LIBS) -o $@
	@echo 'Built: ./calc_qtable.x     (writes ../data/<model>/)'

# Thermal-emission SED of any shipped model, for one
# illuminated cell.  It is built against the T-matrix because the astrodust
# EUV band, when one is asked for, is solved on the oblate spheroid of the Q
# table rather than on a volume-equivalent sphere.
calc_sed.x: $(SRC_PIPE) $(SRC_EUVTM) $(SRC_CLIOPT) $(SRC_CLIAPPLY) src/calc_sed.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $(TMINC) $^ $(TMLIB) $(H5LIBS) -o $@
	@echo 'Built: ./calc_sed.x'

# Size-integrated transport optics (lambda, albedo, <cos>, C_ext/C_abs/C_sca
# per H) for astrodust, DL07, Zubko or any file-defined model -> ../data/.
# Reuses the full SED pipeline, so every table it writes is the same optics
# the RT hosts get from dust_extinction; only the driver differs.
calc_kext.x: $(SRC_PIPE) $(SRC_EUVTM) $(SRC_CLIOPT) src/calc_kext.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $(TMINC) $^ $(TMLIB) $(H5LIBS) -o $@
	@echo 'Built: ./calc_kext.x'

# Does build_dust on the HDF5 products give the same model as the four
# builders on the text products?  Run it after calc_qtable.x and calc_kext.x.
check_build_dust.x: $(SRC_PIPE) $(SRC_EUVTM) src/check_build_dust.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $(TMINC) $^ $(TMLIB) $(H5LIBS) -o $@
	@echo 'Built: ./check_build_dust.x'

# Builds every model from a foreign working directory with an absolute
# data_dir, and checks that a bad one comes back as a status instead of an
# abort.  Links the LIBRARY, not the pipeline sources, because that is what a
# host sees.
check_data_root.x: libsedust.a src/check_data_root.f90
	$(FC) $(FFLAGS) -I. src/check_data_root.f90 libsedust.a $(H5LIBS) -o $@
	@echo 'Built: ./check_data_root.x'

# ---- libsedust.a: the model-agnostic dust-emission library, for a Fortran
# 3D RT code to link.  Usage from the RT:  use dust_lib ;  link
#   -L<sed> -lsedust        with the .mod search path  -I<sed>
# One archive, and no T-matrix in it, so a host that transports no ionizing
# radiation needs nothing from ../tmatrix at all.
#
# A host that DOES want the astrodust EUV band on the spheroid of the Q table
# builds  WITH_TMATRIX=1 make libsedust.a  instead.  The T-matrix objects and
# the spheroid route then go into the same archive (still one -lsedust), and
# the host adds
#       use euv_astrodust_tmatrix, only: use_tmatrix_euv_band_optics
#       call use_tmatrix_euv_band_optics()
# once before build_astrodust.  Without that call, euv_tmatrix = .true. is
# refused with status 6 rather than answered with a sphere.
LIB_SRC = $(SRC_PIPE) src/dust_lib.f90
LIB_OBJ = $(notdir $(LIB_SRC:.f90=.o))

ifdef WITH_TMATRIX

# The T-matrix objects are taken from ../tmatrix's own build, not recompiled
# here: tmatrix_core.f90 needs -fno-inline for its results to be reproducible,
# and that knowledge belongs to the Makefile that owns the source.
libsedust.a: $(LIB_SRC) $(SRC_EUVTM) | $(TMARCH)
	$(FC) $(FFLAGS) $(TMINC) -c $(LIB_SRC) $(SRC_EUVTM)
	rm -rf tmarch.d && mkdir tmarch.d && (cd tmarch.d && ar x ../$(TMARCH))
	rm -f $@
	ar rcs $@ $(LIB_OBJ) $(notdir $(SRC_EUVTM:.f90=.o)) tmarch.d/*.o
	rm -rf tmarch.d
	@echo 'Built: libsedust.a with the T-matrix  (link with -L. -lsedust -I.)'

else

libsedust.a: $(LIB_SRC)
	$(FC) $(FFLAGS) -c $(LIB_SRC)
	rm -f $@
	ar rcs $@ $(LIB_OBJ)
	@echo 'Built: libsedust.a  (link with -L. -lsedust -I.)'

endif

# rt_example: built AS AN EXTERNAL host would -- against libsedust.a and the
# .mod search path, not from the pipeline sources.  Build and run from sed/,
# since the data paths inside are relative to ../.
use_dustlib.x: libsedust.a rt_example/use_dustlib.f90
	$(FC) $(FFLAGS) -I. rt_example/use_dustlib.f90 libsedust.a $(H5LIBS) -o $@
	@echo 'Built: ./use_dustlib.x'

clean:
	rm -f *.o *.mod
	rm -rf tmarch.d

cleanall: clean
	rm -f *.a calc_enthalpy.x calc_sed.x calc_kext.x
	rm -f calc_qtable.x check_build_dust.x use_dustlib.x
	rm -f output/sed_astrodust_*.dat output/enthalpy_*.dat

.PHONY: all clean cleanall

# Every executable recipe compiles the shared module sources into this
# directory, so concurrent targets race on the same .mod files (make -j fails
# with .mod0 rename errors). Serialize until the build uses one object/module
# graph in a separate directory.
.NOTPARALLEL:
