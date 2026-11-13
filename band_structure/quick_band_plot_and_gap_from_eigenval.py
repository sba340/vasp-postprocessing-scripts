#!/usr/bin/env python3
"""Quick band plot and band gap summary from VASP EIGENVAL (non-spin-polarized).

Run this script in a directory that contains EIGENVAL and OUTCAR. It does not
need a KPOINTS file, so the horizontal axis is the k-point index and no
high-symmetry labels are drawn. Use plot_band_structure_from_eigenval.py for a
labelled band structure.

The highest band that has any eigenvalue at or below the Fermi level is taken
as the highest occupied band. The script prints the top of that band, the
bottom of the next band, their difference, the Fermi level, and the middle of
the gap, then saves the plot as a PDF.

Note: the gap is reported as an absolute value, so overlapping bands (a metal)
show a positive number. Check the plot as well as the printed gap.

Usage:
    python quick_band_plot_and_gap_from_eigenval.py -window 4 -outname Band_rough
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument('-window', dest='window', type=float, default=4.0,
                        help='Plot window in eV above and below E_F (default: 4)')
    parser.add_argument('-outname', dest='outname', type=str, default='Band_rough',
                        help='Name of the output PDF without extension')
    return parser.parse_args()


def read_eigenvalues(eigenfile):
    """Return an array of shape (nbands, nkpt)."""
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
        bands.append([float(lines[i].split()[1])
                      for i in range(j, nkpt * block, block)])
    return np.array(bands), nbands, nkpt


def read_fermi_level(outcar):
    """Return the last Fermi level printed in OUTCAR."""
    fermi = None
    with open(outcar) as f:
        for line in f:
            if 'E-fermi' in line:
                fermi = float(line.split()[2])
    if fermi is None:
        sys.exit("No E-fermi line found in OUTCAR.")
    return fermi


def band_edges(eigen, fermi):
    """Return the top of the valence band, bottom of the conduction band,
    the gap, and the middle of the gap."""
    occupied = np.where((eigen <= fermi).any(axis=1))[0]
    if len(occupied) == 0 or occupied.max() + 1 >= eigen.shape[0]:
        sys.exit("Could not find an occupied band with an empty band above it. "
                 "Check that EIGENVAL contains empty bands.")
    homo = occupied.max()
    tvb = np.amax(eigen[homo])
    bcb = np.amin(eigen[homo + 1])
    gap = abs(tvb - bcb)
    middle = tvb + gap / 2.0
    return tvb, bcb, gap, middle


def set_style():
    plt.rc('font', family='serif')
    plt.rc('font', serif=['Times New Roman', 'DejaVu Serif'], size=14)
    plt.rcParams['axes.linewidth'] = 1.20
    for axis in ('x', 'y'):
        plt.rcParams[axis + 'tick.major.size'] = 6
        plt.rcParams[axis + 'tick.minor.size'] = 3
        plt.rcParams[axis + 'tick.major.width'] = 1.20
        plt.rcParams[axis + 'tick.minor.width'] = 1.20
        plt.rcParams[axis + 'tick.direction'] = 'in'


def main():
    args = parse_args()
    for name in ('EIGENVAL', 'OUTCAR'):
        if not os.path.isfile(name):
            sys.exit("%s not found. Required files: EIGENVAL and OUTCAR." % name)

    eigen, nbands, nkpt = read_eigenvalues('EIGENVAL')
    fermi = read_fermi_level('OUTCAR')
    tvb, bcb, gap, middle = band_edges(eigen, fermi)

    print("---------------------------------")
    print("Top of the VB       = ", tvb, "eV")
    print("Bottom of  CB       = ", bcb, "eV")
    print("Band  Gap  is       =  ", round(gap, 4), "eV")
    print("Fermi Energy Ef     = ", fermi)
    print("Band Gap Center BGC = ", middle)
    print("---------------------------------")

    set_style()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for band in eigen:
        ax.plot(band)

    ax.axhline(fermi, color='black', lw=2, ls='--')
    ax.set_xlabel("k-point index", fontsize=16)
    ax.set_ylabel("Energy [eV]", fontsize=16)
    ax.axis([0, nkpt, fermi - args.window, fermi + args.window])
    ax.set_title("Quick band plot")

    plt.savefig(args.outname + '.pdf', bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()

