#!/usr/bin/env python3
"""Plot a band structure and atom-projected band structures with PyProcar.

This script is a command-line wrapper around pyprocar.bandsplot. It produces
one plain band structure and, for each --atoms group given, one band structure
projected onto those atoms.

It uses the PyProcar 5.x interface (bandsplot with the procarfile, outcar and
kpointsfile arguments). PyProcar 6 changed this interface, so with PyProcar 6
install a 5.x release or adapt the call.

Required files (names can be changed with the options below):
    PROCAR   projected band output, from a calculation with LORBIT = 11
    OUTCAR   used for the Fermi level
    KPOINTS  line-mode k-path file

Outputs:
    <prefix>-band.pdf
    <prefix>-band-atoms-<list>.pdf   (one per --atoms group)

Usage:
    python plot_projected_bands_with_pyprocar.py -prefix mysystem -elimit -3 3
    python plot_projected_bands_with_pyprocar.py -prefix mysystem -atoms 1 -atoms 6
    python plot_projected_bands_with_pyprocar.py -atoms 0 1 -atoms 2

Atom numbers follow the PyProcar convention of the installed version, so check
its documentation for whether the first atom is 0 or 1.
"""
import argparse
import inspect
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument('-procar', default='PROCAR', help='PROCAR file (default: PROCAR)')
    parser.add_argument('-outcar', default='OUTCAR', help='OUTCAR file (default: OUTCAR)')
    parser.add_argument('-kpoints', default='KPOINTS', help='KPOINTS file (default: KPOINTS)')
    parser.add_argument('-elimit', nargs=2, type=float, default=[-3.0, 3.0],
                        metavar=('EMIN', 'EMAX'),
                        help='Energy window relative to E_F (default: -3 3)')
    parser.add_argument('-atoms', nargs='+', type=int, action='append', default=[],
                        help='Atom indices for one projected plot; repeat the '
                             'option for several plots')
    parser.add_argument('-color', default='black',
                        help='Line color of the plain band structure (default: black)')
    parser.add_argument('-cmap', default='jet',
                        help='Colormap of the projected plots (default: jet)')
    parser.add_argument('-title', default='Band structure',
                        help='Plot title (default: Band structure)')
    parser.add_argument('-prefix', default='band',
                        help='Prefix of the output file names (default: band)')
    parser.add_argument('-strict', action='store_true',
                        help='Disable permissive PROCAR parsing')
    return parser.parse_args()


def load_pyprocar():
    try:
        import pyprocar
    except ImportError:
        sys.exit("PyProcar is not installed. Install a 5.x release, for example: "
                 "pip install pyprocar==5.6.6")
    if 'procarfile' not in inspect.signature(pyprocar.bandsplot).parameters:
        sys.exit("This script needs the PyProcar 5.x interface of bandsplot "
                 "(procarfile, outcar and kpointsfile arguments). The installed "
                 "version uses a different interface.")
    return pyprocar


def main():
    args = parse_args()
    for path in (args.procar, args.outcar, args.kpoints):
        if not os.path.isfile(path):
            sys.exit("%s not found." % path)
    pyprocar = load_pyprocar()

    common = dict(outcar=args.outcar, kpointsfile=args.kpoints,
                  elimit=args.elimit, grid=True,
                  permissive=not args.strict, show=False)

    pyprocar.bandsplot(args.procar, mode='plain', color=args.color,
                       savefig='%s-band.pdf' % args.prefix,
                       title=args.title, **common)

    for group in args.atoms:
        label = '-'.join(str(a) for a in group)
        pyprocar.bandsplot(args.procar, mode='parametric', cmap=args.cmap,
                           atoms=group,
                           savefig='%s-band-atoms-%s.pdf' % (args.prefix, label),
                           title='%s (atoms %s)' % (args.title, label),
                           **common)


if __name__ == '__main__':
    main()
