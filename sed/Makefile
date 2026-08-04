# SEDust -- SED solver build.
#
# Targets:
#   make                     main_astrodust.x main_dl07.x main_zubko.x
#                            make_enthalpy.x calc_kext.x
#   make libsedust.a         the model-agnostic dust-emission library, with no
#                            T-matrix in it (see below)
#   WITH_TMATRIX=1 make libsedust.a
#                            ... with the T-matrix and the spheroid optics of
#                            the astrodust EUV band in the same archive
#   make calc_kext.x             size-integrated transport optics for any model
#   make calc_polext.x           HD23 polarized extinction per H, vs release
#   make test_scatmat_aligned.x  aligned-scattering reader/interface checks
#   make test_scalar_mode.x      scalar-only build mode checks
#   make test_euv_extension.x    EUV grid-extension (lam_min) checks
#   make use_dustlib_scatmat.x   rt_example aligned-scattering consumer (vs libsedust.a)
#   make clean / cleanall
#
# Run (from this directory; data paths are relative to ../):
#   ./make_enthalpy.x
#   ./main_astrodust.x                 # GD solver (default)
#   ./main_astrodust.x qm              # energy-space, thermal-discrete
#   ./main_astrodust.x qm_dbcon        # energy-space, thermal-continuous
#   ./main_dl07.x
#   ./calc_kext.x astrodust            # -> ../data/kext_astrodust_MW.dat
#   ./calc_kext.x astrodust euv        # carried into the ionizing band
#   ./calc_kext.x dl07
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
# Local variables belong on the run-time stack, not in static storage: the EUV
# band optics of sed_astrodust.f90 are threaded, and they reach
# spheroid_optics.f90 and asymptotic_optics.f90, which this Makefile recompiles
# with its own flags.  The flag is a no-op for GFortran, which already stacks
# local arrays (measured: an object holding a 4096-element local array has an
# empty .bss either way).  It is not a no-op for Intel, where the same object
# carries 32 kB of .bss without it.  Enabling OpenMP happens to imply it on
# both compilers, so nothing is broken today; stating it keeps the reentrancy
# of a threaded path from resting on that side effect.
FFLAGS += -frecursive

# ---- T-matrix library, for the astrodust EUV band only ---------------
# The astrodust model's extreme-ultraviolet band (a lam_min shorter than the Q
# table's own short-wavelength end) is solved on the same oblate spheroid the
# table itself is, which is a T-matrix calculation.  It is the only thing in
# the tree that needs libtmatrix.a, and src/sed_astrodust.f90 reaches it
# through a procedure pointer rather than a use statement, so:
#
#   * `make libsedust.a` links and archives NO T-matrix at all.  A model built
#     from that library answers euv_tmatrix = .true. with status 11.
#   * `WITH_TMATRIX=1 make libsedust.a` adds the T-matrix and the spheroid
#     route to the same archive, so an RT host still links only -lsedust.  The
#     host calls use_tmatrix_euv_band_optics() once before build_astrodust.
#   * the two astrodust CLI drivers (main_astrodust.x, calc_kext.x) always
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
              src/mie.f90 src/q_astrodust.f90 \
              src/q_graphite.f90 src/q_graphite_d16.f90 \
              src/q_graphite_d16_sphere.f90 src/qpah.f90 \
              src/radfield.f90 src/p_sub.f90 \
              src/stoch_qm.f90 $(SRC_DL07MOD) \
              src/scatmat_aligned.f90 src/dust_model_mod.f90 src/zubko_io.f90 \
              src/kext_table.f90 src/sed_astrodust.f90

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

main_astrodust.x: $(SRC_PIPE) $(SRC_EUVTM) src/main_astrodust.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $(TMINC) $^ $(TMLIB) -o $@
	@echo 'Built: ./main_astrodust.x'

# DL07 emission. Its own optics are Mie on the D03 dielectric functions, which
# already run to 6.2e-5 um, so this driver reaches no astrodust EUV band and
# needs no T-matrix.
main_dl07.x: $(SRC_PIPE) src/main_dl07.f90
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./main_dl07.x'

# Zubko/ZDA emission.  Its optics grid (1.0e-3 um, from the DustEM tables)
# reaches far past the illuminated band, so it is the one that exercises the
# difference between "hardest photon the field carries" and "short end of the
# grid".
main_zubko.x: $(SRC_PIPE) src/main_zubko.f90
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./main_zubko.x'

# Size-integrated transport optics (lambda, albedo, <cos>, C_ext/C_abs/C_sca
# per H) for astrodust, DL07, Zubko or any file-defined model -> ../data/.
# Reuses the full SED pipeline and computes the size integral from first
# principles (size_integrated_extinction); the tables it writes here are the
# very ones dust_extinction reads back for an RT host.  The astrodust product
# also carries the dichroic C_polext/H column.
calc_kext.x: $(SRC_PIPE) $(SRC_EUVTM) src/calc_kext.f90 | $(TMARCH)
	$(FC) $(FFLAGS) $(TMINC) $^ $(TMLIB) -o $@
	@echo 'Built: ./calc_kext.x'

# Polarized extinction per H from the orientation-resolved DH21 spheroid
# table, compared against ../data/release/polarized_extinction.dat.
SRC_POLEXT = $(SRC_BASE) src/size_dist.f90 src/q_table_jori.f90 \
             src/calc_polext.f90

calc_polext.x: $(SRC_POLEXT)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./calc_polext.x'

# Aligned-scattering interface checks: reader round-trip, symmetry
# reconstruction, eta algebra, and K consistency against the calc_polext-style
# f_align-weighted size integrals. Reuses the full SED pipeline + dust_lib.
SRC_TSCAT = $(SRC_PIPE) src/dust_lib.f90 src/test_scatmat_aligned.f90

test_scatmat_aligned.x: $(SRC_TSCAT)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./test_scatmat_aligned.x'

# Scalar-only build mode checks: a scalar (load_polarized_optics=.false.) and a
# polarized astrodust build agree on Cext/Cabs/Csca/gbar and the total SED, the
# scalar build carries no polarized optics and opens no polarized Q table, and a
# scalar request combined with an explicit polarized path is a contradiction.
SRC_TSCALAR = $(SRC_PIPE) src/dust_lib.f90 src/test_scalar_mode.f90

test_scalar_mode.x: $(SRC_TSCALAR)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./test_scalar_mode.x'

# EUV grid-extension checks (build_astrodust's optional lam_min): the extended
# grid, the table block left untouched, the join at the Q table's 0.0912 um end,
# and a polarized build whose EUV block reproduces the dichroism and
# birefringence a direct T-matrix calculation of the same spheroid gives.
SRC_TEUV = $(SRC_PIPE) src/dust_lib.f90 src/test_euv_extension.f90

test_euv_extension.x: $(SRC_TEUV)
	$(FC) $(FFLAGS) $^ -o $@
	@echo 'Built: ./test_euv_extension.x'

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
# refused with status 11 rather than answered with a sphere.
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

# rt_example: aligned-scattering Peest-style consumer, built AS AN EXTERNAL host
# would -- against libsedust.a and the .mod search path, not from the pipeline
# sources. Build and run from sed/ (data paths are relative to ../).
use_dustlib_scatmat.x: libsedust.a rt_example/use_dustlib_scatmat.f90
	$(FC) $(FFLAGS) -I. rt_example/use_dustlib_scatmat.f90 libsedust.a -o $@
	@echo 'Built: ./use_dustlib_scatmat.x'

clean:
	rm -f *.o *.mod
	rm -rf tmarch.d

cleanall: clean
	rm -f *.a make_enthalpy.x main_astrodust.x main_dl07.x main_zubko.x \
	      calc_kext.x calc_polext.x \
	      test_scatmat_aligned.x test_scalar_mode.x test_euv_extension.x \
	      use_dustlib_scatmat.x
	rm -f output/astrodust_irem_ours_*.dat output/enthalpy_*.dat

.PHONY: all clean cleanall

# Every executable recipe compiles the shared module sources into this
# directory, so concurrent targets race on the same .mod files (make -j fails
# with .mod0 rename errors). Serialize until the build uses one object/module
# graph in a separate directory.
.NOTPARALLEL:
