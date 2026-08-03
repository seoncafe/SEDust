# SEDust -- reusable T-matrix library and astrodust table driver.
#
# `libtmatrix.a` owns no SEDust dependency and no size-parameter regime
# selection.  `driver/spheroid_optics.f90` owns the regime selection; the
# standalone driver adds only the astrodust dielectric table and grid.

FC = gfortran
BUILD ?= release

# Compiler family.  Only the flag spellings differ; the object list, the
# reentrancy contract, and the tests are the same for both.  On the oneAPI
# releases installed here the `ifort` driver is the IFX compiler, so `ifort`,
# `ifx`, and the MPI wrappers all take the Intel spellings.  Set FC_FAMILY on
# the command line when the compiler is reached through a wrapper whose name
# hides it (mpif90, ftn, ...).
ifneq (,$(filter ifort ifx mpiifort mpiifx,$(notdir $(FC))))
  FC_FAMILY ?= intel
else
  FC_FAMILY ?= gnu
endif

ifeq ($(FC_FAMILY),intel)
  OMPFLAGS ?= -qopenmp
  WARNFLAGS = -warn all
  # -extend-source 132 covers lpd.f, whose only line past column 72 is a
  # comment; keeping it explicit means a future long statement is not silently
  # truncated.  Free-form sources need no length flag because Intel already
  # accepts the whole line.
  LEGACY_BASE = -extend-source 132
  MODERN_BASE = -free
  # -recursive is the Intel spelling of -frecursive: it puts every local
  # variable on the run-time stack, which is what makes a routine reentrant.
  # Intel's default (-auto-scalar) stacks only intrinsic-type scalars and
  # leaves fixed-size local arrays, derived types, and character locals in
  # static storage, so this flag is required, not decorative.
  RECURSIVE_FLAG = -recursive
  # Intel counterpart of -fno-inline.  The reproducibility argument recorded
  # below is a GNU-to-GNU one and does not carry across compilers; the flag is
  # kept so the Intel build has the same "core is not inlined" property.
  NOINLINE_FLAG = -inline-level=0
  # Intel defaults to -fp-model fast, which reassociates and contracts freely.
  # The T-matrix convergence test compares successive iterates against DDELT,
  # so value-unsafe rearrangement moves the iteration count: at the default the
  # 120-case golden regression misses its 1e-10 tolerance and both golden tests
  # fail.  With -fp-model precise every test passes.
  OPT_FLAGS = -O2 -fp-model precise
  LEGACY_QUIET = -warn none
  DEBUG_FLAGS = -O0 -g -check all -traceback -fpe0
else
  OMPFLAGS ?= -fopenmp
  WARNFLAGS = -Wall -Wextra -Wimplicit-interface -Wsurprising
  LEGACY_BASE = -ffixed-line-length-none -fmax-stack-var-size=65536
  MODERN_BASE = -ffree-line-length-none
  RECURSIVE_FLAG = -frecursive
  NOINLINE_FLAG = -fno-inline
  OPT_FLAGS = -O2
  LEGACY_QUIET = -w
  DEBUG_FLAGS = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
endif

ifeq ($(BUILD),debug)
  LEGACY_FFLAGS = $(DEBUG_FLAGS) $(LEGACY_BASE)
  MODERN_FFLAGS = $(DEBUG_FLAGS) $(MODERN_BASE) $(WARNFLAGS)
else
  # Legacy source has pervasive historical diagnostics.  Keep them isolated;
  # modern API/driver/tests always compile with warnings enabled.
  LEGACY_FFLAGS = $(OPT_FLAGS) $(LEGACY_BASE) $(LEGACY_QUIET)
  MODERN_FFLAGS = $(OPT_FLAGS) $(MODERN_BASE) $(WARNFLAGS)
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LDFLAGS = -Wl,-stack_size,0x4000000
else
  LDFLAGS =
endif

PRODUCTION_OBJ = tmatrix_kinds.o tmatrix_status.o tmatrix_types.o tmatrix_core.o \
                 tmatrix_full_direct_bridge.o fixed_orientation_amplitude.o \
                 fixed_orientation_phase_matrix.o scattering_matrix_expansion.o \
                 tmatrix_api.o lpd.o
DRIVER_OBJ = read_index.o asymptotic_optics.o spheroid_optics.o run_tmatrix.o

all: run_tmatrix.x

tmatrix_kinds.o: src/tmatrix_kinds.f90
	$(FC) $(MODERN_FFLAGS) -c -o $@ $<
tmatrix_status.o: src/tmatrix_status.f90
	$(FC) $(MODERN_FFLAGS) -c -o $@ $<
tmatrix_types.o: src/tmatrix_types.f90 tmatrix_kinds.o
	$(FC) $(MODERN_FFLAGS) -I. -c -o $@ $<
tmatrix_core.o: src/tmatrix_core.f90 tmatrix_kinds.o tmatrix_types.o
	# RECURSIVE_FLAG keeps every migrated helper reentrant.  NOINLINE_FLAG is
	# no longer a stack workaround -- TT/TMATR/GSP scratch now lives in the
	# caller's workspace, so an inlined build needs 190 KB of thread stack
	# instead of 3.1 MB and runs safely.  It is kept because inlining changes
	# where the compiler contracts multiply-add pairs, which moves the last
	# bits of the T-matrix convergence test: with gfortran, 2 of 6258 rows of
	# the shipped astrodust table shift by ~1e-6.  Dropping the flag costs
	# about 3% of wall time and that reproducibility; both are deliberate.
	# Reproducibility is only claimed within one compiler.
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) $(NOINLINE_FLAG) -I. -c -o $@ $<
tmatrix_full_direct_bridge.o: src/tmatrix_full_direct_bridge.f90 tmatrix_kinds.o tmatrix_types.o tmatrix_core.o
	# Consumers may invoke the production API from a threading runtime even
	# when they do not compile with OpenMP.  Keep this bridge reentrant without
	# introducing an OpenMP link dependency.
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
fixed_orientation_amplitude.o: src/fixed_orientation_amplitude.f90 tmatrix_kinds.o tmatrix_types.o
	# AMPL's CAL(NPN4,NPN4) and the four VIGAMPL work vectors are automatic
	# locals.  RECURSIVE_FLAG keeps them on the caller's stack so concurrent
	# callers with distinct workspaces do not share them.
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
fixed_orientation_phase_matrix.o: src/fixed_orientation_phase_matrix.f90 tmatrix_kinds.o
	# Consumers may build the phase matrix from a threading runtime; keep it
	# reentrant like the amplitude routine that feeds it.
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
scattering_matrix_expansion.o: src/scattering_matrix_expansion.f90 tmatrix_kinds.o tmatrix_types.o
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
tmatrix_api.o: src/tmatrix_api.f90 tmatrix_kinds.o tmatrix_status.o tmatrix_types.o tmatrix_core.o tmatrix_full_direct_bridge.o fixed_orientation_amplitude.o fixed_orientation_phase_matrix.o
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
lpd.o: src/lpd.f
	$(FC) $(LEGACY_FFLAGS) $(RECURSIVE_FLAG) -c -o $@ $<

libtmatrix.a: $(PRODUCTION_OBJ)
	ar rcs $@ $^
	@echo 'Built production archive: ./libtmatrix.a'

# driver layer.  The reentrancy contract covers it too: the SED solver calls
# spheroid_q from an OpenMP loop, and spheroid_q calls rayleigh_limit and
# geometric_optics_limit, so the locals of all three must live on the calling
# thread's stack.  Under gfortran that already happens by default; Intel's
# default (-auto-scalar) stacks only intrinsic-type scalars and sends fixed-size
# local arrays, derived types, and character locals to static storage, where
# concurrent callers overwrite each other.  Measured with nm on ifx 2025.3
# without the flag, asymptotic_optics.o puts the dS/nhat/shat/phat vectors of
# fresnel_opaque_absorption and the qe_o/qa_o arrays of geometric_optics_limit
# in .bss, and test_parallel_spheroid_optics fails from 2 threads up; with the
# flag the object has no .bss.  Which locals a given release chooses to stack is
# not something to bet the contract on, so the flag goes on every driver object,
# not only on src/.
DRIVER_FFLAGS = $(MODERN_FFLAGS) $(RECURSIVE_FLAG)

read_index.o: driver/read_index.f90
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
asymptotic_optics.o: driver/asymptotic_optics.f90 tmatrix_kinds.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
spheroid_optics.o: driver/spheroid_optics.f90 tmatrix_kinds.o tmatrix_status.o tmatrix_api.o asymptotic_optics.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
run_tmatrix.o: driver/run_tmatrix.f90 tmatrix_api.o tmatrix_status.o read_index.o asymptotic_optics.o spheroid_optics.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
run_tmatrix.x: libtmatrix.a $(DRIVER_OBJ)
	$(FC) $(DRIVER_OBJ) -L. -ltmatrix $(LDFLAGS) -o $@
	@echo 'Built: ./run_tmatrix.x'
	@echo 'Smoke test:  TMATRIX_OUTPUT_FILE=/tmp/table.dat ./run_tmatrix.x test'
	@echo 'Full sweep:  ./run_tmatrix.x'

evaluate_one.x: libtmatrix.a examples/evaluate_one.f90
	$(FC) $(MODERN_FFLAGS) -I. examples/evaluate_one.f90 -L. -ltmatrix $(LDFLAGS) -o $@
example: evaluate_one.x

test_cases.o: tests/test_cases.f90
	# This immutable case catalogue is read concurrently by the OpenMP-only
	# harnesses; make its module procedures reentrant without linking OpenMP.
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -c -o $@ $<
test_api_regression.x: libtmatrix.a test_cases.o tests/test_api_regression.f90
	$(FC) $(MODERN_FFLAGS) -I. tests/test_api_regression.f90 test_cases.o -L. -ltmatrix $(LDFLAGS) -o $@
test_api_status.x: libtmatrix.a tests/test_api_status.f90
	$(FC) $(MODERN_FFLAGS) -I. tests/test_api_status.f90 -L. -ltmatrix $(LDFLAGS) -o $@
test_full_direct_regression.x: libtmatrix.a test_cases.o tests/test_full_direct_regression.f90
	$(FC) $(MODERN_FFLAGS) -I. tests/test_full_direct_regression.f90 test_cases.o -L. -ltmatrix $(LDFLAGS) -o $@
test_parallel_api.x: libtmatrix.a test_cases.o tests/test_parallel_api.f90
	$(FC) $(MODERN_FFLAGS) $(OMPFLAGS) -I. tests/test_parallel_api.f90 test_cases.o -L. -ltmatrix $(LDFLAGS) -o $@
test_parallel_full_direct.x: libtmatrix.a test_cases.o tests/test_parallel_full_direct.f90
	$(FC) $(MODERN_FFLAGS) $(OMPFLAGS) -I. tests/test_parallel_full_direct.f90 test_cases.o -L. -ltmatrix $(LDFLAGS) -o $@
test_parallel_scattering_matrix.x: libtmatrix.a test_cases.o tests/test_parallel_scattering_matrix.f90
	$(FC) $(MODERN_FFLAGS) $(OMPFLAGS) -I. tests/test_parallel_scattering_matrix.f90 test_cases.o -L. -ltmatrix $(LDFLAGS) -o $@
# The only harness that enters the driver layer: it drives spheroid_q through
# all three size-parameter regimes, so the Rayleigh and geometric-optics paths
# through asymptotic_optics are covered as well.  It carries RECURSIVE_FLAG
# itself because it evaluates one case in an internal procedure whose locals
# are not covered by an OpenMP private clause; the harness must be reentrant
# for the objects under test to be the only thing it measures.
test_parallel_spheroid_optics.x: libtmatrix.a asymptotic_optics.o spheroid_optics.o tests/test_parallel_spheroid_optics.f90
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) $(OMPFLAGS) -I. tests/test_parallel_spheroid_optics.f90 asymptotic_optics.o spheroid_optics.o -L. -ltmatrix $(LDFLAGS) -o $@

test: test-production-package test-default-backend test_api_regression.x test_api_status.x test_full_direct_regression.x test_parallel_api.x test_parallel_full_direct.x test_parallel_scattering_matrix.x test_parallel_spheroid_optics.x
	./test_api_regression.x
	./test_api_status.x
	./test_full_direct_regression.x
	OMP_NUM_THREADS=4 ./test_parallel_api.x
	OMP_NUM_THREADS=4 ./test_parallel_full_direct.x
	OMP_NUM_THREADS=4 ./test_parallel_scattering_matrix.x
	OMP_NUM_THREADS=4 ./test_parallel_spheroid_optics.x

test-default-backend:
	@awk 'index($$0, "subroutine tmatrix_eval("), index($$0, "end subroutine tmatrix_eval")' src/tmatrix_api.f90 | grep -q 'call full_direct_eval'
	@if awk 'index($$0, "subroutine tmatrix_eval("), index($$0, "end subroutine tmatrix_eval")' src/tmatrix_api.f90 | grep -Eq 'legacy_eval|ctctt_eval|cbess_eval|!\$$omp[[:space:]]+critical'; then \
		echo 'default tmatrix_eval must not enter a legacy or serialized backend'; exit 1; \
	fi
	@echo 'PASS: default tmatrix_eval selects only the full-direct backend'

test-production-package: libtmatrix.a example
	@if nm -u libtmatrix.a | grep -Eq 'GOMP_|__kmpc_|omp_'; then \
		echo 'production archive must not require an OpenMP runtime'; exit 1; \
	fi
	@if ar t libtmatrix.a | grep -Eq '(^|/)tmd_one|tmatrix_(legacy|ctctt|cbess)_bridge'; then \
		echo 'production archive must not contain legacy or staged bridge objects'; exit 1; \
	fi
	@echo 'PASS: production archive is OpenMP- and legacy-bridge-free'

test-full-direct-parallel: test_parallel_full_direct.x
	OMP_NUM_THREADS=1 ./test_parallel_full_direct.x
	OMP_NUM_THREADS=2 ./test_parallel_full_direct.x
	OMP_NUM_THREADS=4 ./test_parallel_full_direct.x
	OMP_NUM_THREADS=8 ./test_parallel_full_direct.x

test-api-parallel: test_parallel_api.x
	OMP_NUM_THREADS=1 ./test_parallel_api.x
	OMP_NUM_THREADS=2 ./test_parallel_api.x
	OMP_NUM_THREADS=4 ./test_parallel_api.x
	OMP_NUM_THREADS=8 ./test_parallel_api.x

test-scattering-matrix-parallel: test_parallel_scattering_matrix.x
	OMP_NUM_THREADS=1 ./test_parallel_scattering_matrix.x
	OMP_NUM_THREADS=2 ./test_parallel_scattering_matrix.x
	OMP_NUM_THREADS=4 ./test_parallel_scattering_matrix.x
	OMP_NUM_THREADS=8 ./test_parallel_scattering_matrix.x

test-spheroid-optics-parallel: test_parallel_spheroid_optics.x
	OMP_NUM_THREADS=1 ./test_parallel_spheroid_optics.x
	OMP_NUM_THREADS=2 ./test_parallel_spheroid_optics.x
	OMP_NUM_THREADS=4 ./test_parallel_spheroid_optics.x
	OMP_NUM_THREADS=8 ./test_parallel_spheroid_optics.x

# The driver normally writes output/.  Direct it to a caller-supplied temporary
# path when smoke testing so generated artifacts never enter the repository.
smoke: run_tmatrix.x
	TMATRIX_OUTPUT_FILE=/tmp/tmatrix_smoke_table.dat ./run_tmatrix.x test
	cmp /tmp/tmatrix_smoke_table.dat tests/golden/standalone_test.dat
	rm -f /tmp/tmatrix_smoke_table.dat

test-debug:
	$(MAKE) clean
	$(MAKE) BUILD=debug test

clean:
	rm -f *.o *.mod *.a *.x

cleanall: clean
	rm -f run_tmatrix.x evaluate_one.x test_api_regression.x test_api_status.x test_full_direct_regression.x test_parallel_api.x test_parallel_full_direct.x test_parallel_scattering_matrix.x test_parallel_spheroid_optics.x

.PHONY: all clean cleanall example test smoke test-debug test-default-backend test-production-package test-full-direct-parallel test-api-parallel test-scattering-matrix-parallel test-spheroid-optics-parallel

# Module files are emitted into this directory.  Serialize until the project
# has a dedicated object/module build tree.
.NOTPARALLEL:
