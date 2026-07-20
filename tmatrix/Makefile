# SEDust -- T-matrix optical properties of randomly oriented oblate spheroids.
#
# Two executables share the same T-matrix core:
#   run_tmatrix.x   cross-section table (Q_ext, Q_abs, Q_sca, albedo, g)
#   run_scatmat.x   size-integrated scattering (Mueller) matrix
#
# Common sources:
#   src/tmd_one.f                 TMD_ONE / TMD_ONE_SCATMAT wrappers plus the
#                                 Mishchenko inner subroutines
#   src/lpd.f                     LAPACK routines
#   ../sed/src/constants.f90      working precision
#   driver/read_index.f90         dielectric-function reader
#   driver/asymptotic_optics.f90  Rayleigh and geometric-optics limits
#
# Run:
#   make                     # build both executables
#   ./run_tmatrix.x test     # 49-point sanity sweep (~seconds)
#   ./run_tmatrix.x          # full 169 x 1129 = 191k sweep (1-10 hours)
#   ./run_scatmat.x test     # 3-wavelength scattering-matrix smoke test
#   ./run_scatmat.x 0.55     # scattering matrix at one wavelength
#   make clean
#
# The full sweep writes output/q_astrodust_P0.20_Fe0.00_1.400.dat, which the
# SED solver reads.  That table ships with SEDust, so regenerating it is
# optional -- and `make cleanall` deliberately does NOT delete it.

FC = gfortran

# Notes on flags:
#  -fmax-stack-var-size=65536: arrays > 64 KB go to heap, not stack.  Without
#    this the Mishchenko COMMON storage (>16 MB) can overflow the stack.
FFLAGS = -O2 -ffixed-line-length-none -ffree-line-length-none -w \
         -fmax-stack-var-size=65536
# debug build: replace with -O0 -g -fcheck=all -fbacktrace -Wall

# macOS needs an explicit 64 MB main-thread stack: its 8 MB default is too
# small for Mishchenko's local arrays and ulimit -s cannot raise it.  GNU ld
# on Linux neither accepts -stack_size nor needs it.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LDFLAGS = -Wl,-stack_size,0x4000000
else
  LDFLAGS =
endif

SRC_TM   = src/tmd_one.f src/lpd.f src/ampl_oriented.f
SRC_COMM = ../sed/src/constants.f90 \
           driver/read_index.f90 driver/asymptotic_optics.f90
SRC_Q    = $(SRC_COMM) driver/run_tmatrix.f90
SRC_QJ   = $(SRC_COMM) driver/tmatrix_oriented.f90 driver/run_q_jori.f90
SRC_SM   = $(SRC_COMM) ../sed/src/size_dist.f90 driver/run_scatmat.f90
SRC_CMP  = $(SRC_COMM) driver/compare_rayleigh_jori.f90
SRC_CMP2 = $(SRC_COMM) driver/tmatrix_oriented.f90 \
           driver/compare_tmatrix_jori.f90
SRC_CMP3 = $(SRC_COMM) driver/tmatrix_oriented.f90 \
           driver/compare_tmatrix_abs.f90
SRC_CMP4 = $(SRC_COMM) driver/tmatrix_oriented.f90 \
           driver/compare_go_jori.f90
SRC_CMP5 = $(SRC_COMM) driver/tmatrix_oriented.f90 \
           driver/compare_birefringence.f90

all: run_tmatrix.x run_scatmat.x run_q_jori.x

# Each executable is compiled in its own directory so that the two builds do
# not fight over .mod files of the shared modules.
run_tmatrix.x: $(SRC_TM) $(SRC_Q)
	@mkdir -p build/q
	$(FC) $(FFLAGS) -Isrc -Jbuild/q $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./run_tmatrix.x'
	@echo 'Smoke test:  ./run_tmatrix.x test'
	@echo 'Full sweep:  ./run_tmatrix.x'

run_scatmat.x: $(SRC_TM) $(SRC_SM)
	@mkdir -p build/sm
	$(FC) $(FFLAGS) -Isrc -Jbuild/sm $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./run_scatmat.x'
	@echo 'Smoke test:  ./run_scatmat.x test'

# Orientation-resolved Q table for the polarization channel: sweeps the DH21
# grid through oriented_cross_sections and writes the file q_table_jori.f90
# reads.  Shares the T-matrix core and the asymptotic limits with the two
# executables above.
run_q_jori.x: $(SRC_TM) $(SRC_QJ)
	@mkdir -p build/qj
	$(FC) $(FFLAGS) -Isrc -Jbuild/qj $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./run_q_jori.x'
	@echo 'Smoke test:  ./run_q_jori.x test'
	@echo 'Full sweep:  ./run_q_jori.x'

# Stage-A check: oriented Rayleigh Q_abs vs the HD23 orientation-resolved
# table over the small-x subset.  No T-matrix core is needed.
compare_rayleigh_jori.x: $(SRC_CMP)
	@mkdir -p build/cmp
	$(FC) $(FFLAGS) -Isrc -Jbuild/cmp $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_rayleigh_jori.x'
	@echo 'Run:   ./compare_rayleigh_jori.x'

# Stage-B check: oriented T-matrix Q_ext vs the HD23 orientation-resolved
# table over the 0.1 < x < 50 subset.  Uses the T-matrix core plus the
# fixed-orientation amplitude routine AMPL (src/ampl_oriented.f).
compare_tmatrix_jori.x: $(SRC_TM) $(SRC_CMP2)
	@mkdir -p build/cmp2
	$(FC) $(FFLAGS) -Isrc -Jbuild/cmp2 $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_tmatrix_jori.x'
	@echo 'Run:   ./compare_tmatrix_jori.x'

# Stage-C check: oriented T-matrix Q_abs (= Q_ext - Q_sca, with Q_sca from
# the fixed-orientation phase-matrix integral) vs the HD23 orientation-
# resolved table over the 0.1 < x < 50 subset.  Uses the T-matrix core plus
# the fixed-orientation amplitude routine AMPL (src/ampl_oriented.f).
compare_tmatrix_abs.x: $(SRC_TM) $(SRC_CMP3)
	@mkdir -p build/cmp3
	$(FC) $(FFLAGS) -Isrc -Jbuild/cmp3 $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_tmatrix_abs.x'
	@echo 'Run:   ./compare_tmatrix_abs.x'

# Geometric-optics check: oriented Q_ext (projected-area extinction paradox)
# and Q_abs (opaque-grain Fresnel surface integral) vs the HD23 orientation-
# resolved table over the x > 50 subset, plus the continuity checks across
# x = 50 that use the T-matrix core (src/ampl_oriented.f).
compare_go_jori.x: $(SRC_TM) $(SRC_CMP4)
	@mkdir -p build/cmp4
	$(FC) $(FFLAGS) -Isrc -Jbuild/cmp4 $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_go_jori.x'
	@echo 'Run:   ./compare_go_jori.x'

# Birefringence certification: the real part of the forward-amplitude
# difference (circular-polarization optic) reconstructed from the T-matrix,
# checked internally against the imaginary part by Kramers-Kronig, plus the
# Rayleigh-analytic, continuity, and stored-precision anchors.
compare_birefringence.x: $(SRC_TM) $(SRC_CMP5)
	@mkdir -p build/cmp5
	$(FC) $(FFLAGS) -Isrc -Jbuild/cmp5 $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_birefringence.x'
	@echo 'Run:   ./compare_birefringence.x'

clean:
	rm -f *.o *.mod
	rm -rf build

cleanall: clean
	rm -f run_tmatrix.x run_scatmat.x run_q_jori.x compare_rayleigh_jori.x \
	      compare_tmatrix_jori.x compare_tmatrix_abs.x compare_go_jori.x \
	      compare_birefringence.x

.PHONY: all clean cleanall
