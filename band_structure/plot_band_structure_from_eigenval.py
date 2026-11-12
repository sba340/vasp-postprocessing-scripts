#!/usr/bin/env python3
"""Plot a band structure from VASP output (non-spin-polarized calculation).

Run this script in a directory that contains:
    EIGENVAL  eigenvalues along the k-path
    KPOINTS   line-mode k-path file
    OUTCAR    used to read the Fermi level and the k-point coordinates

The KPOINTS file is expected in line-mode format, with the high-symmetry
label as the last item of each coordinate line, for example:

    0.0 0.0 0.0 ! G
    0.5 0.0 0.0 ! X

Outputs:
    <outname>.pdf   band structure plotted relative to the Fermi level
    BAND.dat        k-path distance and eigenvalues, one block per band

Usage:
    python plot_band_structure_from_eigenval.py -emin -3 -emax 3 -outname Band_VASP
"""
import argparse
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument('-emin', dest='emin', type=float, default=-3.0,
                        help='Lower energy limit relative to E_F (default: -3)')
    parser.add_argument('-emax', dest='emax', type=float, default=3.0,
                        help='Upper energy limit relative to E_F (default: 3)')
    parser.add_argument('-outname', dest='outname', type=str, default='Band_VASP',
                        help='Name of the output file without extension')
    return parser.parse_args()


def check_inputs(files=('OUTCAR', 'KPOINTS', 'EIGENVAL')):
    for name in files:
        if not os.path.isfile(name):
            sys.exit("%s not found. Required files: EIGENVAL, KPOINTS "
                     "(line-mode path) and OUTCAR." % name)


def read_eigenvalues(eigenfile):
    """Return a list of bands, each a list of eigenvalues along the k-path."""
    with open(eigenfile) as f:
        for _ in range(5):
            f.readline()
        header = f.readline().split()
        nkpt = int(header[1])
        nbands = int(header[2])
        f.readline()
        f.readline()
        lines = f.readlines()

    block = nbands + 2  # band lines plus a blank line and the next k-point line
    bands = []
    for j in range(nbands):
        band = [float(lines[i].split()[1])
                for i in range(j, nkpt * block, block)]
        bands.append(band)
    return bands, nbands


def read_kpath_labels(kpointfile):
    """Return the high-symmetry labels from a line-mode KPOINTS file."""
    with open(kpointfile) as f:
        f.readline()
        points_per_segment = int(f.readline().split()[0])
        f.readline()
        f.readline()
        contents = f.read().rstrip().split('\n')

    names = [contents[i].split()[-1].strip('!#')
             for i in range(0, len(contents), 3)]
    last_label = contents[-1].split()[-1].strip('!#')
    names.append(last_label)
    names = [r'$\Gamma$' if n in ('G', r'\Gamma') else n for n in names]
    return names, points_per_segment


def read_outcar(outcar):
    """Return the Fermi level and the k-point coordinates from OUTCAR."""
    with open(outcar) as f:
        lines = f.readlines()

    fermi = None
    for line in lines:
        if 'E-fermi' in line:
            fermi = float(line.split()[2])

    started = "k-points in units of 2pi/SCALE"
    ended = "k-points in reciprocal lattice and weights"
    kpoints = []
    reading = False
    for line in lines:
        if re.search(started, line):
            reading = True
            continue
        if re.search(ended, line):
            break
        if reading and line.strip():
            kpoints.append([float(x) for x in line.split()[:3]])
    return fermi, np.array(kpoints)


def k_distances(kpoints, points_per_segment):
    """Cumulative distance along the path and at each high-symmetry point."""
    kdist = [0.0]
    total = 0.0
    for i in range(len(kpoints) - 1):
        total += np.linalg.norm(kpoints[i + 1] - kpoints[i])
        kdist.append(total)
    kdist_hsym = [kdist[i] for i in range(0, len(kdist), points_per_segment)]
    kdist_hsym.append(kdist[-1])
    return kdist, kdist_hsym


def set_style():
    plt.rc('font', family='serif')
    plt.rc('font', serif=['Times New Roman', 'DejaVu Serif'], size=20)
    plt.rcParams['axes.linewidth'] = 1.20
    for axis in ('x', 'y'):
        plt.rcParams[axis + 'tick.major.size'] = 6
        plt.rcParams[axis + 'tick.minor.size'] = 3
        plt.rcParams[axis + 'tick.major.width'] = 1.20
        plt.rcParams[axis + 'tick.minor.width'] = 1.20
        plt.rcParams[axis + 'tick.direction'] = 'in'


def main():
    args = parse_args()
    check_inputs()

    bands, nbands = read_eigenvalues('EIGENVAL')
    kname, points_per_segment = read_kpath_labels('KPOINTS')
    fermi, kpoints = read_outcar('OUTCAR')
    kdist, kdist_hsym = k_distances(kpoints, points_per_segment)

    set_style()
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    with open('BAND.dat', 'w') as fout:
        for band in bands:
            energies = np.array(band)
            ax.plot(kdist, energies - fermi, lw=1.5, c='navy')
            for i in range(len(kdist)):
                fout.write("%8d %18.8f %18.8f\n" % (i + 1, kdist[i], energies[i]))
            fout.write("\n")

    ax.axhline(0, color='red', lw=1, ls='-.')
    ax.set_ylabel(r"$E-E_{\mathrm{F}}$ [eV]", fontsize=18)
    ax.axis([0, max(kdist), args.emin, args.emax])

    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9, pad=0.2)
    ax.annotate(args.outname, xy=(0.02, args.emax * 0.82), xycoords='data',
                horizontalalignment='left', verticalalignment='bottom',
                bbox=bbox_props)

    for x in kdist_hsym:
        ax.axvline(x, color='gray', lw=1, ls='--')

    plt.xticks(kdist_hsym, kname)
    plt.savefig(args.outname + '.pdf', bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


if __name__ == '__main__':
    main()

