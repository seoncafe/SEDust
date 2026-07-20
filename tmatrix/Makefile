# SEDust -- T-matrix cross sections for randomly oriented oblate spheroids.
#
# Builds run_tmatrix.x by linking:
#   src/tmd_one.f            TMD_ONE wrapper + Mishchenko inner subroutines
#   src/lpd.f                LAPACK routines
#   ../sed/src/constants.f90 working precision
#   driver/read_index.f90    dielectric-function reader
#   driver/asymptotic_optics.f90  analytic Rayleigh / geometric-optics limits
#   driver/run_tmatrix.f90   driver
#
# Run:
#   make                     # build run_tmatrix.x
#   ./run_tmatrix.x test     # 49-point sanity sweep (~seconds)
#   ./run_tmatrix.x          # full 169 x 1129 = 191k sweep (1-10 hours)
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

SRC_TM  = src/tmd_one.f src/lpd.f
SRC_F90 = ../sed/src/constants.f90 \
          driver/read_index.f90 driver/asymptotic_optics.f90 driver/run_tmatrix.f90

run_tmatrix.x: $(SRC_TM) $(SRC_F90)
	$(FC) $(FFLAGS) -Isrc $^ $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./run_tmatrix.x'
	@echo 'Smoke test:  ./run_tmatrix.x test'
	@echo 'Full sweep:  ./run_tmatrix.x'

clean:
	rm -f *.o *.mod

cleanall: clean
	rm -f run_tmatrix.x

.PHONY: clean cleanall
