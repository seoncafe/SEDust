#!/usr/bin/env python3
"""Store the regime codes of an existing SEDust HDF5 product as one byte each.

WHAT THE POLICY IS.  A product holds two kinds of dataset, and the storage type
follows what the numbers ARE (the same rule the Fortran writers in
sed/src/sedust_h5.f90 apply when a product is generated fresh):

  values  Q_abs, Q_sca, Q_ext, Q_re, g, Z, F_tot, F_ref, the C_* cross
          sections, K_abs, albedo, and the coordinate axes lambda, a_eff,
          theta, theta_i, theta_s, phi, band_lambda.  All 64-bit float, so a
          number read back out of a product is the number that was computed.
  codes   regime.  Three distinct values -- 0 T-matrix, 10 Rayleigh dipole,
          20 geometric optics -- so one signed byte each.  Eight bytes to
          carry three states says the quantity is continuous when it is not,
          and h5_read_2d_int reads either storage back as integers.

Only the astrodust product has a regime dataset; every other product is
already on the policy and is reported as such without being touched.

WHY IN-PLACE CONVERSION AND NOT REGENERATION.  The products hold results that
cost hours to recompute -- the aligned scattering matrix /polarized/Z of the
astrodust product alone is 81 MB of T-matrix output -- and rerunning
calc_qtable.x REPLACES the file, taking /polarized with it.  Nothing about the
storage type of a code requires recomputing anything, so this rewrites what is
on disk instead.

WHY REWRITE AND NOT MODIFY.  HDF5 never reclaims the space of an unlinked or
replaced dataset; it leaves a hole.  A new file written dataset by dataset
avoids that.  Everything else is carried over unchanged: shape, chunking,
compression, shuffle, and every dataset, group and file attribute.  Every
dataset other than regime is required to come out bit-for-bit identical, and
the run fails if one does not.

Usage:
    python3 pyutil/migrate_h5_regime_int8.py [--dry-run] [FILE ...]

With no FILE given it converts every data/<model>/sedust_<model>.h5 under the
tree this script lives in.  Each file that needs converting is written to a
temporary neighbour, verified against the original, and only then moved into
place.
"""

import argparse
import os
import sys

import h5py
import numpy as np

# The dataset of codes.  Matched on the dataset's base name, as the products
# name it under the group of the component it belongs to.
CODE_NAMES = frozenset(('regime',))

CODE_DTYPE = np.dtype('int8')


def is_code(name):
    return name.rsplit('/', 1)[-1] in CODE_NAMES


def target_dtype(name, dset):
    if is_code(name):
        return CODE_DTYPE
    return dset.dtype


def code_datasets(path):
    """(every code dataset in the file, those not already one byte)."""
    found, pending = [], []

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset) and is_code(name):
            found.append(name)
            if obj.dtype != CODE_DTYPE:
                pending.append(name)

    with h5py.File(path, 'r') as f:
        f.visititems(visit)
    return found, pending


def copy_attrs(src, dst):
    for key, val in src.attrs.items():
        dst.attrs.create(key, val)


def convert(src_path, dst_path):
    """Write src_path to dst_path with the code datasets stored as one byte."""
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

        n_code = n_val = 0
        for name in sorted(set(names_a) & set(names_b)):
            da, db = names_a[name], names_b[name]
            want = target_dtype(name, da)
            if db.dtype != want:
                ok = False
                lines.append('    %s: dtype %s, expected %s' % (name, db.dtype, want))
            if da.attrs.keys() != db.attrs.keys():
                ok = False
                lines.append('    %s: attributes changed' % name)

            va, vb = da[...], db[...]
            if va.shape != vb.shape:
                ok = False
                lines.append('    %s: shape %s -> %s' % (name, va.shape, vb.shape))
                continue

            if is_code(name):
                n_code += 1
                # The codes have to survive exactly, and they have to have been
                # whole numbers in the first place -- a non-integral entry
                # would mean this dataset is not codes at all.
                fa = np.asarray(va, dtype=np.float64)
                if not np.all(fa == np.rint(fa)):
                    ok = False
                    lines.append('    %s: not integral, so not a code dataset' % name)
                if not np.array_equal(np.asarray(vb, dtype=np.float64), fa):
                    ok = False
                    lines.append('    %s: code values changed' % name)
                lines.append('    %s: %s -> %s, distinct codes %s'
                             % (name, da.dtype, db.dtype,
                                sorted(int(x) for x in np.unique(fa))))
            else:
                n_val += 1
                # Everything that is not a code must come through untouched:
                # compare the raw bytes, not the values, so that a -0.0 or a
                # NaN payload that changed would still be caught.
                if va.dtype != vb.dtype or va.tobytes() != vb.tobytes():
                    ok = False
                    lines.append('    %s: NOT BIT-IDENTICAL' % name)

        lines.append('    values bit-identical: %d;  code datasets converted: %d'
                     % (n_val, n_code))

        if dict(a['/'].attrs) != dict(b['/'].attrs):
            ok = False
            lines.append('    file attributes changed')
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
        print('%s' % path)
        found, pending = code_datasets(path)
        if not found:
            print('    no code datasets; nothing to convert')
            continue
        if not pending:
            print('    codes already one byte; not rewritten')
            continue
        before = os.path.getsize(path)
        tmp = path + '.i8.tmp'
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
            print('    converted')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
