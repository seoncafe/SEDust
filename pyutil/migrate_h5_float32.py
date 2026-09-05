#!/usr/bin/env python3
"""Bring an existing SEDust HDF5 product onto the 32-bit storage policy.

WHAT THE POLICY IS.  Three classes of dataset, decided by what the numbers ARE
rather than by how large the array is (the same rule the Fortran writers in
sed/src/sedust_h5.f90 apply when a product is generated fresh):

  quantities  Q_abs, Q_sca, Q_ext, Q_re, g, Z, F_tot, F_ref, the C_* cross
              sections, K_abs, albedo -- anything that came out of a
              calculation.  Stored as 32-bit.  The physics carries its own
              uncertainty far above float32's 1.2e-7 relative resolution: the
              dielectric functions are good to three figures, the 1/3-2/3
              graphite orientation average is 40% wrong at 10 um, and the
              radiative transfer that consumes these has Monte Carlo noise of
              1e-2 to 1e-3.
  axes        lambda, a_eff, theta, theta_i, theta_s, phi, band_lambda.  Left
              at 64 bits.  Their problem is not accuracy but DISTINCTNESS and
              monotonicity -- the astrodust wavelength grid resolves each
              X-ray absorption edge with a pair of points 6.7e-7 apart in
              relative wavelength, only six to eleven representable float32
              values apart -- and they are 0.06% of the payload.
  codes       regime.  One signed byte, already.

WHY IN-PLACE CONVERSION AND NOT REGENERATION.  The products hold results that
cost hours to recompute -- the aligned scattering matrix /polarized/Z of the
astrodust product alone is 81 MB of T-matrix output -- and rerunning
calc_qtable.x REPLACES the file, taking /polarized with it.  Nothing about the
storage type of a number requires recomputing the number, so this rewrites what
is on disk instead.

WHY REWRITE AND NOT MODIFY.  HDF5 never reclaims the space of an unlinked or
replaced dataset; it leaves a hole.  A new file written dataset by dataset is
what actually makes the file smaller.  Everything else is carried over
unchanged: shape, chunking, compression, shuffle, and every dataset, group and
file attribute.

Usage:
    python3 pyutil/migrate_h5_float32.py [--dry-run] [FILE ...]

With no FILE given it migrates every data/<model>/sedust_<model>.h5 under the
tree this script lives in.  Each file is converted to a temporary neighbour,
verified against the original (dataset set unchanged, storage classes correct,
values agreeing to within float32 resolution), and only then moved into place.
"""

import argparse
import os
import sys

import h5py
import numpy as np

# Coordinate axes: the datasets that stay 64-bit.  Matched on the dataset's
# base name, since the same axis appears under several groups.
AXIS_NAMES = frozenset((
    'lambda', 'a_eff', 'theta', 'theta_i', 'theta_s', 'phi', 'band_lambda',
))

# float32 has a 24-bit significand, so the worst relative error of a correctly
# rounded conversion is 2**-24 = 5.96e-8.  Allow a shade over twice that.
F32_EPS = 2.0 ** -24
TOL = 1.3e-7

# The smallest NORMAL float32.  Below it the format keeps the value but loses
# significand bits, and below about 1.4e-45 it has nothing left and the value
# becomes zero.  Two of the products reach down there: the longest-wavelength
# end of kext/C_sca is 6.2e-40 cm^2/H, and the aligned phase matrix Z carries
# entries as small as 4.5e-65 um^2/sr/H against a peak of 3.9e-14 -- fifty
# orders of magnitude below the peak of the very array they sit in, in a table
# that is only ever summed or interpolated linearly, never divided by.
#
# So agreement is judged two ways, and both must hold:
#   rel_normal  elementwise relative error over the elements at or above
#               FLT_MIN_NORMAL, where float32 is fully relatively accurate.
#               This is the float32-resolution test proper.
#   rel_peak    the largest absolute change divided by the largest magnitude in
#               the dataset.  This is what says an underflowed entry did not
#               matter: a number 1e-30 of the peak contributes nothing to any
#               integral over the array.
FLT_MIN_NORMAL = 1.1754943508222875e-38


def storage_class(name, dset):
    """'axis', 'code' or 'quantity' -- what this dataset is."""
    if dset.dtype.kind in 'iu':
        return 'code'
    if name.rsplit('/', 1)[-1] in AXIS_NAMES:
        return 'axis'
    return 'quantity'


def target_dtype(name, dset):
    cls = storage_class(name, dset)
    if cls == 'quantity' and dset.dtype == np.float64:
        return np.dtype('<f4')
    return dset.dtype


def copy_attrs(src, dst):
    for key, val in src.attrs.items():
        dst.attrs.create(key, val)


def convert(src_path, dst_path):
    """Write src_path to dst_path with the quantity datasets narrowed."""
    with h5py.File(src_path, 'r') as src, h5py.File(dst_path, 'w') as dst:
        copy_attrs(src, dst)

        def visit(name, obj):
            if isinstance(obj, h5py.Group):
                copy_attrs(obj, dst.require_group(name))
            elif isinstance(obj, h5py.Dataset):
                out = dst.create_dataset(
                    name,
                    shape=obj.shape,
                    dtype=target_dtype(name, obj),
                    chunks=obj.chunks,
                    compression=obj.compression,
                    compression_opts=obj.compression_opts,
                    shuffle=obj.shuffle,
                    fletcher32=obj.fletcher32,
                )
                out[...] = obj[...]
                copy_attrs(obj, out)
            else:
                raise TypeError('unhandled object %r at %s' % (type(obj), name))

        src.visititems(visit)


def verify(src_path, dst_path):
    """Return (ok, list of report lines)."""
    lines = []
    ok = True
    with h5py.File(src_path, 'r') as a, h5py.File(dst_path, 'r') as b:
        names_a, names_b = {}, {}
        a.visititems(lambda n, o: names_a.__setitem__(n, o)
                     if isinstance(o, h5py.Dataset) else None)
        b.visititems(lambda n, o: names_b.__setitem__(n, o)
                     if isinstance(o, h5py.Dataset) else None)

        if set(names_a) != set(names_b):
            ok = False
            lines.append('    DATASET SET CHANGED: only in original %s; only in new %s'
                         % (sorted(set(names_a) - set(names_b)),
                            sorted(set(names_b) - set(names_a))))
        lines.append('    datasets: %d -> %d' % (len(names_a), len(names_b)))

        worst, worst_name = 0.0, ''
        n_axis = n_code = n_q = 0
        subnormal = []
        for name in sorted(set(names_a) & set(names_b)):
            da, db = names_a[name], names_b[name]
            cls = storage_class(name, da)
            want = target_dtype(name, da)
            if db.dtype != want:
                ok = False
                lines.append('    %s: dtype %s, expected %s' % (name, db.dtype, want))
            if cls == 'axis':
                n_axis += 1
                if db.dtype != np.float64:
                    ok = False
                    lines.append('    AXIS %s is not float64 (%s)' % (name, db.dtype))
            elif cls == 'code':
                n_code += 1
            else:
                n_q += 1
            if da.attrs.keys() != db.attrs.keys():
                ok = False
                lines.append('    %s: attributes changed' % name)

            va = np.asarray(da[...], dtype=np.float64)
            vb = np.asarray(db[...], dtype=np.float64)
            if va.shape != vb.shape:
                ok = False
                lines.append('    %s: shape %s -> %s' % (name, va.shape, vb.shape))
                continue
            if va.size == 0:
                continue
            diff = np.abs(vb - va)
            peak = np.max(np.abs(va))
            normal = np.abs(va) >= FLT_MIN_NORMAL
            rel_normal = (np.max(diff[normal] / np.abs(va[normal]))
                          if np.any(normal) else 0.0)
            rel_peak = diff.max() / peak if peak > 0.0 else 0.0
            n_sub = int(np.sum((va != 0.0) & ~normal))
            if n_sub:
                subnormal.append((name, n_sub, va.size, peak))
            rel = max(rel_normal, rel_peak)
            if rel > worst:
                worst, worst_name = rel, name
            if rel_normal > TOL:
                ok = False
                lines.append('    %s: relative change %.3e over the normal range '
                             'exceeds %.1e' % (name, rel_normal, TOL))
            if rel_peak > TOL:
                ok = False
                lines.append('    %s: change %.3e of the dataset peak exceeds %.1e'
                             % (name, rel_peak, TOL))

        lines.append('    axes kept float64: %d;  code datasets: %d;  '
                     'quantities now float32: %d' % (n_axis, n_code, n_q))
        lines.append('    worst relative change: %.3e  (%s)   [float32 eps = %.2e]'
                     % (worst, worst_name or '-', F32_EPS))
        for name, n_sub, n_tot, peak in subnormal:
            lines.append('    note: %s has %d of %d entries below float32\'s normal '
                         'range (dataset peak %.3e); they keep at most a few '
                         'significant bits, and the smallest become zero'
                         % (name, n_sub, n_tot, peak))

        for grp_name in ('/',):
            if dict(a[grp_name].attrs) != dict(b[grp_name].attrs):
                ok = False
                lines.append('    file attributes changed')
        # group attributes
        groups_a, groups_b = {}, {}
        a.visititems(lambda n, o: groups_a.__setitem__(n, dict(o.attrs))
                     if isinstance(o, h5py.Group) else None)
        b.visititems(lambda n, o: groups_b.__setitem__(n, dict(o.attrs))
                     if isinstance(o, h5py.Group) else None)
        if set(groups_a) != set(groups_b):
            ok = False
            lines.append('    group set changed')
        for g in sorted(set(groups_a) & set(groups_b)):
            if list(groups_a[g]) != list(groups_b[g]):
                ok = False
                lines.append('    group %s: attributes changed' % g)
    return ok, lines


def default_products(root):
    data = os.path.join(root, 'data')
    found = []
    if os.path.isdir(data):
        for model in sorted(os.listdir(data)):
            path = os.path.join(data, model, 'sedust_%s.h5' % model)
            if os.path.isfile(path):
                found.append(path)
    return found


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('files', nargs='*', metavar='FILE')
    ap.add_argument('--dry-run', action='store_true',
                    help='convert and verify, but leave the original in place')
    args = ap.parse_args(argv)

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    paths = args.files or default_products(root)
    if not paths:
        print('no products found under %s/data' % root)
        return 1

    failed = 0
    for path in paths:
        before = os.path.getsize(path)
        tmp = path + '.f32.tmp'
        print('%s' % path)
        convert(path, tmp)
        ok, lines = verify(path, tmp)
        after = os.path.getsize(tmp)
        for line in lines:
            print(line)
        print('    size: %.2f MB -> %.2f MB  (%.1f%%)'
              % (before / 1e6, after / 1e6, 100.0 * after / before))
        if not ok:
            failed += 1
            print('    FAILED -- original left untouched, new file at %s' % tmp)
        elif args.dry_run:
            print('    dry run -- new file left at %s' % tmp)
        else:
            os.replace(tmp, path)
            print('    migrated')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
