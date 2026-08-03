# Upstream fixed-form sources

This directory preserves the original Mishchenko standalone source snapshots:

- `tmd.lp.f` and `tmd.par.f` implement the historical polydisperse,
  randomly-oriented T-matrix program.
- `ampld.lp.f` and `ampld.par.f` implement the historical amplitude and phase
  matrix program.

They are provenance material, are never compiled, and are not part of
`libtmatrix.a` or the standalone SEDust table driver build. `tmd.par.f` is the
historical array-limit include; `src/tmatrix_core.f90` owns the same limits as
private parameters.
