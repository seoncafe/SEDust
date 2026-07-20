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

SRC_TM   = src/tmd_one.f src/lpd.f
SRC_COMM = ../sed/src/constants.f90 \
           driver/read_index.f90 driver/asymptotic_optics.f90
SRC_Q    = $(SRC_COMM) driver/run_tmatrix.f90
SRC_SM   = $(SRC_COMM) ../sed/src/size_dist.f90 driver/run_scatmat.f90
SRC_CMP  = $(SRC_COMM) driver/compare_rayleigh_jori.f90

all: run_tmatrix.x run_scatmat.x

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

# Stage-A check: oriented Rayleigh Q_abs vs the HD23 orientation-resolved
# table over the small-x subset.  No T-matrix core is needed.
compare_rayleigh_jori.x: $(SRC_CMP)
	@mkdir -p build/cmp
	$(FC) $(FFLAGS) -Isrc -Jbuild/cmp $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_rayleigh_jori.x'
	@echo 'Run:   ./compare_rayleigh_jori.x'

clean:
	rm -f *.o *.mod
	rm -rf build

cleanall: clean
	rm -f run_tmatrix.x run_scatmat.x compare_rayleigh_jori.x

.PHONY: all clean cleanall
