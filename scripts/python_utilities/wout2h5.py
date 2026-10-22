#!/usr/bin/env python3
"""Convert a VMEC wout file (netCDF) into the HDF5 file read by GYACOMO's vmec module.

Usage: python wout2h5.py wout_xxx.nc [wout_xxx.h5] [--gzip N]

Only the arrays needed by src/vmec_mod.f90 are copied, with VMEC's own names.
Scalars are written as 1-element datasets so that futils' getarr can read them.
2D Fourier arrays keep the netCDF (radius, mn_mode) layout, which Fortran sees
as (mnmax, ns), the layout used by LIBSTELL's read_wout_mod and by stella.

wout files are classic netCDF-3 so scipy.io.netcdf_file reads them without a
netCDF library; netCDF4 is used as a fallback for netCDF-4/HDF5 based files.
"""
import argparse
import os
import sys

import numpy as np
import h5py

SCALARS_INT = ['nfp', 'ns', 'mnmax', 'mnmax_nyq', 'mpol', 'ntor', 'signgs']
SCALARS_REAL = ['Aminor_p', 'Rmajor_p']
ARRAYS_1D = ['xm', 'xn', 'xm_nyq', 'xn_nyq', 'phi', 'phips', 'iotas', 'iotaf', 'presf']
ARRAYS_2D = ['rmnc', 'zmns', 'lmns', 'bmnc', 'gmnc',
             'bsupumnc', 'bsupvmnc', 'bsubumnc', 'bsubvmnc', 'bsubsmns']
ARRAYS_2D_ASYM = ['rmns', 'zmnc', 'lmnc', 'bmns', 'gmns',
                  'bsupumns', 'bsupvmns', 'bsubumns', 'bsubvmns', 'bsubsmnc']


def open_wout(path):
    """Return a dict-like {name: numpy array} for the wout file."""
    with open(path, 'rb') as fh:
        magic = fh.read(3)
    if magic == b'CDF':
        from scipy.io import netcdf_file
        nc = netcdf_file(path, 'r', mmap=False)
        return {k: np.array(v.data) for k, v in nc.variables.items()}
    import netCDF4  # noqa: E402 (optional dependency)
    nc = netCDF4.Dataset(path, 'r')
    return {k: np.array(v[:]) for k, v in nc.variables.items()}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('wout', help='VMEC wout_*.nc file')
    ap.add_argument('out', nargs='?', help='output HDF5 file (default: same name, .h5)')
    ap.add_argument('--gzip', type=int, default=4, help='gzip level for 2D arrays (0 = none)')
    args = ap.parse_args(argv)

    out = args.out or os.path.splitext(args.wout)[0] + '.h5'
    v = open_wout(args.wout)

    lasym = int(v.get('lasym__logical__', v.get('lasym', 0)))
    missing = [n for n in SCALARS_INT + SCALARS_REAL + ARRAYS_1D + ARRAYS_2D if n not in v]
    if lasym:
        missing += [n for n in ARRAYS_2D_ASYM if n not in v]
    if missing:
        sys.exit('wout2h5: missing variables in %s: %s' % (args.wout, ', '.join(missing)))

    kw = dict(compression='gzip', compression_opts=args.gzip) if args.gzip > 0 else {}
    with h5py.File(out, 'w') as h:
        h.attrs['source_file'] = os.path.basename(args.wout)
        h.attrs['wout2h5_version'] = 1
        for n in SCALARS_INT:
            h.create_dataset(n, data=np.array([int(v[n])], dtype=np.int32))
        h.create_dataset('lasym', data=np.array([lasym], dtype=np.int32))
        for n in SCALARS_REAL:
            h.create_dataset(n, data=np.array([float(v[n])], dtype=np.float64))
        for n in ARRAYS_1D:
            h.create_dataset(n, data=np.asarray(v[n], dtype=np.float64))
        for n in ARRAYS_2D + (ARRAYS_2D_ASYM if lasym else []):
            h.create_dataset(n, data=np.asarray(v[n], dtype=np.float64), **kw)
    print('wout2h5: wrote %s (nfp=%d, ns=%d, mnmax=%d, mnmax_nyq=%d, lasym=%d)'
          % (out, int(v['nfp']), int(v['ns']), int(v['mnmax']), int(v['mnmax_nyq']), lasym))


if __name__ == '__main__':
    main()
