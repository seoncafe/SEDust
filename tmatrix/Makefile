# SEDust -- T-matrix optical properties of astrodust spheroids, random and
# fixed orientation.
#
# `libtmatrix.a` is the reusable T-matrix library: caller-owned workspace, no
# COMMON, no SEDust dependency, and no size-parameter regime selection.  Its
# sources in src/ are kept byte-identical to SEDust_v1.00/tmatrix/src/, and the
# unmodified Mishchenko standalone sources they were migrated from are kept,
# uncompiled, in reference/upstream/.
#
# driver/ owns everything the library deliberately does not:
#   spheroid_optics.f90         size-parameter regime selection (Rayleigh
#                               dipole / T-matrix / geometric optics), the IERR
#                               redirection, and the flag convention
#   asymptotic_optics.f90       the two closed-form limits, orientation-resolved
#   tmatrix_oriented.f90        orientation-resolved Q_ext/Q_sca/Q_abs from the
#                               forward amplitude and the scattering-sphere
#                               integral
#   scattering_matrix_oriented.f90  fixed-orientation Mueller matrix
#   aligned_population_optics.f90   size-integrated aligned-population products
#   read_index.f90              dielectric-function reader
#
# Executables:
#   run_tmatrix.x            random-orientation cross-section table
#   run_scatmat.x            size-integrated random-orientation Mueller matrix
#   run_q_jori.x             orientation-resolved Q table (polarization channel)
#   run_scatmat_aligned.x    size-integrated aligned Mueller matrix (OpenMP)
#   compare_*.x              verification drivers
#   euv_polarized_optics.x   orientation-resolved optics shortward of 0.0912 um
#
# Run:
#   make                     # run_tmatrix.x run_scatmat.x run_q_jori.x
#   make all-compare         # every verification driver
#   ./run_tmatrix.x test     # 49-point sanity sweep (~seconds)
#   ./run_tmatrix.x          # full 169 x 1129 = 191k sweep (1-10 hours)
#   ./run_scatmat.x test     # 3-wavelength scattering-matrix smoke test
#   make clean
#
# The full sweeps write output/q_astrodust_*.dat, which the SED solver reads.
# Those tables ship with SEDust, so regenerating them is optional -- and
# `make cleanall` deliberately does NOT delete them.

FC = gfortran
BUILD ?= release

# Compiler family.  Only the flag spellings differ; the object list and the
# reentrancy contract are the same for both.  On the oneAPI releases installed
# here the `ifort` driver is the IFX compiler, so `ifort`, `ifx`, and the MPI
# wrappers all take the Intel spellings.  Set FC_FAMILY on the command line
# when the compiler is reached through a wrapper whose name hides it (mpif90,
# ftn, ...).
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
  # Intel counterpart of -fno-inline, kept so the Intel build has the same
  # "core is not inlined" property.  The reproducibility argument recorded at
  # the tmatrix_core.o rule is a GNU-to-GNU one and does not carry across
  # compilers.
  NOINLINE_FLAG = -inline-level=0
  # Intel defaults to -fp-model fast, which reassociates and contracts freely.
  # The T-matrix convergence test compares successive iterates against DDELT,
  # so value-unsafe rearrangement moves the iteration count.  -fp-model precise
  # is what keeps the solver on the iteration path the golden data records.
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
  # The one remaining fixed-form source (LAPACK) carries pervasive historical
  # diagnostics.  Keep them isolated; library, driver, and programs always
  # compile with warnings enabled.
  LEGACY_FFLAGS = $(OPT_FLAGS) $(LEGACY_BASE) $(LEGACY_QUIET)
  MODERN_FFLAGS = $(OPT_FLAGS) $(MODERN_BASE) $(WARNFLAGS)
endif

# macOS needs an explicit 64 MB main-thread stack: its 8 MB default is too
# small for the local arrays and ulimit -s cannot raise it.  GNU ld on Linux
# neither accepts -stack_size nor needs it.
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LDFLAGS = -Wl,-stack_size,0x4000000
else
  LDFLAGS =
endif

# ---------------------------------------------------------------------------
# libtmatrix.a -- identical sources and identical flags to SEDust_v1.00.
# ---------------------------------------------------------------------------
PRODUCTION_OBJ = tmatrix_kinds.o tmatrix_status.o tmatrix_types.o tmatrix_core.o \
                 tmatrix_full_direct_bridge.o fixed_orientation_amplitude.o \
                 fixed_orientation_phase_matrix.o scattering_matrix_expansion.o \
                 tmatrix_api.o lpd.o

tmatrix_kinds.o: src/tmatrix_kinds.f90
	$(FC) $(MODERN_FFLAGS) -c -o $@ $<
tmatrix_status.o: src/tmatrix_status.f90
	$(FC) $(MODERN_FFLAGS) -c -o $@ $<
tmatrix_types.o: src/tmatrix_types.f90 tmatrix_kinds.o
	$(FC) $(MODERN_FFLAGS) -I. -c -o $@ $<
tmatrix_core.o: src/tmatrix_core.f90 tmatrix_kinds.o tmatrix_types.o
	# RECURSIVE_FLAG keeps every migrated helper reentrant.  NOINLINE_FLAG is
	# kept because inlining changes where the compiler contracts multiply-add
	# pairs, which moves the last bits of the T-matrix convergence test.
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) $(NOINLINE_FLAG) -I. -c -o $@ $<
tmatrix_full_direct_bridge.o: src/tmatrix_full_direct_bridge.f90 tmatrix_kinds.o tmatrix_types.o tmatrix_core.o
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
fixed_orientation_amplitude.o: src/fixed_orientation_amplitude.f90 tmatrix_kinds.o tmatrix_types.o
	# AMPL's CAL(NPN4,NPN4) and the four VIGAMPL work vectors are automatic
	# locals.  RECURSIVE_FLAG keeps them on the caller's stack so concurrent
	# callers sharing one workspace do not share them.
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
fixed_orientation_phase_matrix.o: src/fixed_orientation_phase_matrix.f90 tmatrix_kinds.o
	# Called from the OpenMP node loop of aligned_population_optics, so keep
	# it reentrant like the amplitude routine that feeds it.
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
scattering_matrix_expansion.o: src/scattering_matrix_expansion.f90 tmatrix_kinds.o tmatrix_types.o
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
tmatrix_api.o: src/tmatrix_api.f90 tmatrix_kinds.o tmatrix_status.o tmatrix_types.o tmatrix_core.o tmatrix_full_direct_bridge.o fixed_orientation_amplitude.o fixed_orientation_phase_matrix.o
	$(FC) $(MODERN_FFLAGS) $(RECURSIVE_FLAG) -I. -c -o $@ $<
lpd.o: src/lpd.f
	$(FC) $(LEGACY_FFLAGS) $(RECURSIVE_FLAG) -c -o $@ $<

libtmatrix.a: $(PRODUCTION_OBJ)
	ar rcs $@ $^
	@echo 'Built library archive: ./libtmatrix.a'

# ---------------------------------------------------------------------------
# driver layer.  Everything carries RECURSIVE_FLAG: the Mueller-matrix and
# dipole routines are called from the OpenMP node loop of
# aligned_population_optics, and the SED solver calls spheroid_q -- and through
# it rayleigh_limit and geometric_optics_limit -- from an OpenMP loop of its
# own, so their locals must live on the calling thread's stack.  Under gfortran
# that already happens by default; Intel's default (-auto-scalar) stacks only
# intrinsic-type scalars and leaves local arrays, derived types, and character
# locals in static storage, where concurrent callers overwrite each other.
# ---------------------------------------------------------------------------
DRIVER_FFLAGS = $(MODERN_FFLAGS) $(RECURSIVE_FLAG)

read_index.o: driver/read_index.f90
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
asymptotic_optics.o: driver/asymptotic_optics.f90 tmatrix_kinds.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
spheroid_optics.o: driver/spheroid_optics.f90 tmatrix_kinds.o tmatrix_status.o tmatrix_api.o asymptotic_optics.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
tmatrix_oriented.o: driver/tmatrix_oriented.f90 tmatrix_api.o tmatrix_status.o asymptotic_optics.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
scattering_matrix_oriented.o: driver/scattering_matrix_oriented.f90 tmatrix_api.o asymptotic_optics.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
aligned_population_optics.o: driver/aligned_population_optics.f90 tmatrix_api.o tmatrix_status.o asymptotic_optics.o scattering_matrix_oriented.o
	# The Z-node loop of one size is the only parallel region in this tree.
	$(FC) $(DRIVER_FFLAGS) $(OMPFLAGS) -I. -c -o $@ $<

# Readers shared with the SED solver.  size_dist.f90 fixes its own working
# precision and pulls in nothing else; q_table_jori.f90 takes the HD23
# alignment law from it.
size_dist.o: ../sed/src/size_dist.f90
	$(FC) $(DRIVER_FFLAGS) -c -o $@ $<
q_table_jori.o: ../sed/src/q_table_jori.f90 size_dist.o
	$(FC) $(DRIVER_FFLAGS) -c -o $@ $<
# ../sed/src/mie.f90 takes its working precision from ../sed/src/constants.f90.
# That module is compiled for euv_polarized_optics.x alone, which is the only
# target needing the Mie comparison; nothing else in this directory reads it.
constants.o: ../sed/src/constants.f90
	$(FC) $(DRIVER_FFLAGS) -c -o $@ $<
mie.o: ../sed/src/mie.f90 constants.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<

# ---------------------------------------------------------------------------
# executables
# ---------------------------------------------------------------------------
OBJ_BASE     = read_index.o asymptotic_optics.o
OBJ_ORIENTED = $(OBJ_BASE) tmatrix_oriented.o
OBJ_ALIGNED  = $(OBJ_ORIENTED) scattering_matrix_oriented.o size_dist.o \
               q_table_jori.o aligned_population_optics.o

all: run_tmatrix.x run_scatmat.x run_q_jori.x

all-compare: compare_rayleigh_jori.x compare_tmatrix_jori.x compare_tmatrix_abs.x \
             compare_go_jori.x compare_birefringence.x compare_scatmat_aligned.x

run_tmatrix.o: driver/run_tmatrix.f90 tmatrix_api.o tmatrix_status.o read_index.o asymptotic_optics.o spheroid_optics.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
run_tmatrix.x: libtmatrix.a $(OBJ_BASE) spheroid_optics.o run_tmatrix.o
	$(FC) $(OBJ_BASE) spheroid_optics.o run_tmatrix.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./run_tmatrix.x'
	@echo 'Smoke test:  ./run_tmatrix.x test'
	@echo 'Full sweep:  ./run_tmatrix.x   then   make lyman_cut'

# The two scalar Q tables that ship in output/ are one calculation.  The sweep
# above covers the whole DH21_wave_to_12keV axis and writes the EUV product;
# this target writes the other one, by dropping every wavelength shortward of
# the Lyman limit from it.  Row selection only -- no value is recomputed -- so
# the two agree cell for cell wherever both carry a row, and the wavelength
# count in the header is recounted rather than assumed.  (The
# orientation-resolved table is not part of this: it stays on the
# 1129-wavelength DH21_wave axis and has nothing to cut.)
Q_EUV   = output/q_astrodust_P0.20_Fe0.00_1.400_euv.dat
Q_LYMAN = output/q_astrodust_P0.20_Fe0.00_1.400.dat
lyman_cut: $(Q_EUV)
	@n=`awk '$$0 !~ /^#/ && NF >= 8 && $$1+0 >= 0.0912 { print $$1 }' $(Q_EUV) \
	     | sort -u -g | wc -l`; \
	awk -v n=$$n 'NR==1 { print; next } \
	     /^#/ { if ($$0 ~ /lambda stride/) sub(/of [0-9]+/, "of " n); \
	            if ($$0 ~ /^#   lambda\[um\]/) { \
	              print "# Wavelengths shortward of 0.0912 um (the Lyman limit) dropped from"; \
	              print "#   q_astrodust_P0.20_Fe0.00_1.400_euv.dat.  Row selection only: no value is"; \
	              print "#   recomputed, so every row here is the row that file carries." } \
	            print; next } \
	     (NF >= 8 && $$1+0 >= 0.0912) { print }' $(Q_EUV) > $(Q_LYMAN); \
	echo "Wrote $(Q_LYMAN)  ($$n wavelengths, 0.0912 um and longward)"

run_scatmat.o: driver/run_scatmat.f90 tmatrix_api.o tmatrix_status.o scattering_matrix_expansion.o read_index.o asymptotic_optics.o size_dist.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
run_scatmat.x: libtmatrix.a $(OBJ_BASE) size_dist.o run_scatmat.o
	$(FC) $(OBJ_BASE) size_dist.o run_scatmat.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./run_scatmat.x'
	@echo 'Smoke test:  ./run_scatmat.x test'

# Orientation-resolved Q table for the polarization channel: sweeps the DH21
# grid through oriented_cross_sections and writes the file q_table_jori.f90
# reads.
run_q_jori.o: driver/run_q_jori.f90 tmatrix_api.o tmatrix_status.o read_index.o tmatrix_oriented.o asymptotic_optics.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
run_q_jori.x: libtmatrix.a $(OBJ_ORIENTED) run_q_jori.o
	$(FC) $(OBJ_ORIENTED) run_q_jori.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./run_q_jori.x'
	@echo 'Smoke test:  ./run_q_jori.x test'
	@echo 'Full sweep:  ./run_q_jori.x'

# Size-integrated fixed-orientation (aligned) scattering matrix.  The Z-node
# loop is parallelized with OpenMP; the threads share one workspace, which AMPL
# only reads.
run_scatmat_aligned.o: driver/run_scatmat_aligned.f90 tmatrix_api.o tmatrix_status.o scattering_matrix_expansion.o read_index.o size_dist.o q_table_jori.o aligned_population_optics.o
	$(FC) $(DRIVER_FFLAGS) $(OMPFLAGS) -I. -c -o $@ $<
run_scatmat_aligned.x: libtmatrix.a $(OBJ_ALIGNED) run_scatmat_aligned.o
	$(FC) $(OMPFLAGS) $(OBJ_ALIGNED) run_scatmat_aligned.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./run_scatmat_aligned.x'
	@echo 'Smoke test:  ./run_scatmat_aligned.x test'

# Stage-A check: oriented Rayleigh Q_abs vs the HD23 orientation-resolved
# table over the small-x subset.  No T-matrix solve is needed.
compare_rayleigh_jori.o: driver/compare_rayleigh_jori.f90 tmatrix_kinds.o read_index.o asymptotic_optics.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
compare_rayleigh_jori.x: libtmatrix.a $(OBJ_BASE) compare_rayleigh_jori.o
	$(FC) $(OBJ_BASE) compare_rayleigh_jori.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_rayleigh_jori.x'

# Stage-B check: oriented T-matrix Q_ext vs the HD23 orientation-resolved
# table over the 0.1 < x < 50 subset.
compare_tmatrix_jori.o: driver/compare_tmatrix_jori.f90 tmatrix_api.o tmatrix_status.o read_index.o tmatrix_oriented.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
compare_tmatrix_jori.x: libtmatrix.a $(OBJ_ORIENTED) compare_tmatrix_jori.o
	$(FC) $(OBJ_ORIENTED) compare_tmatrix_jori.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_tmatrix_jori.x'

# Stage-C check: oriented T-matrix Q_abs (= Q_ext - Q_sca, with Q_sca from the
# fixed-orientation phase-matrix integral) vs the HD23 table.
compare_tmatrix_abs.o: driver/compare_tmatrix_abs.f90 tmatrix_api.o tmatrix_status.o read_index.o asymptotic_optics.o tmatrix_oriented.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
compare_tmatrix_abs.x: libtmatrix.a $(OBJ_ORIENTED) compare_tmatrix_abs.o
	$(FC) $(OBJ_ORIENTED) compare_tmatrix_abs.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_tmatrix_abs.x'

# Geometric-optics check: oriented Q_ext (projected-area extinction paradox)
# and Q_abs (opaque-grain Fresnel surface integral) vs the HD23 table over the
# x > 50 subset, plus the continuity checks across x = 50.
compare_go_jori.o: driver/compare_go_jori.f90 tmatrix_api.o tmatrix_status.o read_index.o asymptotic_optics.o tmatrix_oriented.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
compare_go_jori.x: libtmatrix.a $(OBJ_ORIENTED) compare_go_jori.o
	$(FC) $(OBJ_ORIENTED) compare_go_jori.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_go_jori.x'

# Birefringence certification: the real part of the forward-amplitude
# difference, checked against its imaginary part by Kramers-Kronig, plus the
# Rayleigh-analytic, continuity, and stored-precision anchors.
compare_birefringence.o: driver/compare_birefringence.f90 tmatrix_api.o tmatrix_status.o read_index.o asymptotic_optics.o tmatrix_oriented.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
compare_birefringence.x: libtmatrix.a $(OBJ_ORIENTED) compare_birefringence.o
	$(FC) $(OBJ_ORIENTED) compare_birefringence.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_birefringence.x'

# Fixed-orientation scattering-matrix engine and aligned population: closure,
# optical theorem, random-orientation average, Rayleigh limit, symmetry (A-E),
# size-integrated closure vs the jori table (F), the Rayleigh K sin^2 law (G),
# and the meridional rotations at general geometry (H).  Anchor F calls the
# size-integrated accumulation, so this links OpenMP like run_scatmat_aligned.x.
compare_scatmat_aligned.o: driver/compare_scatmat_aligned.f90 tmatrix_api.o tmatrix_status.o scattering_matrix_expansion.o read_index.o tmatrix_oriented.o asymptotic_optics.o scattering_matrix_oriented.o size_dist.o q_table_jori.o aligned_population_optics.o
	$(FC) $(DRIVER_FFLAGS) $(OMPFLAGS) -I. -c -o $@ $<
compare_scatmat_aligned.x: libtmatrix.a $(OBJ_ALIGNED) compare_scatmat_aligned.o
	$(FC) $(OMPFLAGS) $(OBJ_ALIGNED) compare_scatmat_aligned.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./compare_scatmat_aligned.x'

# Orientation-resolved optics shortward of the Q table's 0.0912 um end, at
# arbitrary (lambda, a_eff) nodes: the Rayleigh dipole limit, the T-matrix with
# its status IERR reported per node, and the geometric-optics limit side by
# side, plus the volume-equivalent Mie sphere.
euv_polarized_optics.o: driver/euv_polarized_optics.f90 tmatrix_api.o tmatrix_status.o read_index.o asymptotic_optics.o tmatrix_oriented.o mie.o
	$(FC) $(DRIVER_FFLAGS) -I. -c -o $@ $<
euv_polarized_optics.x: libtmatrix.a $(OBJ_ORIENTED) constants.o mie.o euv_polarized_optics.o
	$(FC) $(OBJ_ORIENTED) constants.o mie.o euv_polarized_optics.o -L. -ltmatrix $(LDFLAGS) -o $@
	@echo
	@echo 'Built: ./euv_polarized_optics.x'
	@echo 'Run:   ./euv_polarized_optics.x NODEFILE OUTFILE [sca|nosca] [XTM_HI]'

# The smallest possible external consumer of the public API: no SEDust module,
# no input file, no workspace beyond the one it initializes itself.  It is what
# a host code links against, so it is built the same way a host would.
evaluate_one.x: libtmatrix.a examples/evaluate_one.f90
	$(FC) $(MODERN_FFLAGS) -I. examples/evaluate_one.f90 -L. -ltmatrix $(LDFLAGS) -o $@
example: evaluate_one.x

clean:
	rm -f *.o *.mod *.a
	rm -rf build

cleanall: clean
	rm -f run_tmatrix.x run_scatmat.x run_q_jori.x compare_rayleigh_jori.x \
	      compare_tmatrix_jori.x compare_tmatrix_abs.x compare_go_jori.x \
	      compare_birefringence.x compare_scatmat_aligned.x run_scatmat_aligned.x \
	      euv_polarized_optics.x evaluate_one.x

.PHONY: all all-compare example clean cleanall lyman_cut

# Module files are emitted into this directory.  Serialize until the project
# has a dedicated object/module build tree.
.NOTPARALLEL:
