#!/usr/bin/env python3
"""Plot a band structure next to the total density of states from VASP output.

Run this script in a directory that contains:
    EIGENVAL  eigenvalues along the k-path
    KPOINTS   line-mode k-path file
    OUTCAR    used to read the Fermi level and the k-point coordinates
    DOSCAR    total density of states

Spin polarization is detected from the number of columns in the total DOS
block of DOSCAR (3 columns: non-spin-polarized; 5 columns: spin-polarized).
For spin-polarized runs the spin-down bands are drawn in a separate panel on
the left, the spin-up bands in the middle, and the DOS on the right with the
spin-down DOS drawn as negative values. EIGENVAL is read with the same
setting, so both the older (2 or 3 column) and the newer (3 or 5 column)
EIGENVAL layouts are handled.

The Fermi level used to align the bands and the DOS can be taken from OUTCAR
or from DOSCAR (-fermi). With the default (auto), spin-polarized runs use the
DOSCAR value and non-spin-polarized runs use the OUTCAR value. Use the same
setting for runs that you want to compare.

Output:
    <outname>.pdf  (default: Band_DOS_VASP.pdf)

Usage:
    python plot_band_structure_with_dos.py -emin -3 -emax 3 -dosmax 10
"""
import argparse
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument('-emin', dest='emin', type=float, default=None,
                        help='Lower energy limit relative to E_F '
                             '(default: -4 spin-polarized, -3 otherwise)')
    parser.add_argument('-emax', dest='emax', type=float, default=None,
                        help='Upper energy limit relative to E_F '
                             '(default: 4 spin-polarized, 3 otherwise)')
    parser.add_argument('-dosmax', dest='dosmax', type=float, default=10.0,
                        help='Half-width of the DOS axis in states/eV '
                             '(default: 10)')
    parser.add_argument('-fermi', dest='fermi', choices=['auto', 'outcar', 'doscar'],
                        default='auto',
                        help='Source of the Fermi level (default: auto)')
    parser.add_argument('-title', dest='title', type=str,
                        default='Band structure and DOS',
                        help='Figure title')
    parser.add_argument('-outname', dest='outname', type=str,
                        default='Band_DOS_VASP',
                        help='Output file name without extension')
    return parser.parse_args()


def check_inputs(files=('OUTCAR', 'DOSCAR', 'KPOINTS', 'EIGENVAL')):
    for name in files:
        if not os.path.isfile(name):
            sys.exit("%s not found. Required files: EIGENVAL, KPOINTS "
                     "(line-mode path), OUTCAR (Fermi level and k-points) and "
                     "DOSCAR (density of states)." % name)


def close_polygon(values):
    """Pad with zeros so that filled areas close at the baseline."""
    return np.insert(np.append(values, 0), 0, 0)


def read_total_dos(doscar):
    """Return (energy, dos_up, dos_down, e_fermi, is_magnetic).

    dos_down is None for a non-spin-polarized calculation.
    """
    with open(doscar) as f:
        for _ in range(5):
            f.readline()
        header = f.readline().split()
        no_steps = int(header[2])
        e_fermi = float(header[3])
        dos = [[float(x) for x in f.readline().split()] for _ in range(no_steps)]

    ncols = len(dos[0])
    if ncols == 5:
        is_magnetic = True
        print("DOSCAR is magnetic")
    elif ncols == 3:
        is_magnetic = False
        print("DOSCAR is non-magnetic")
    else:
        sys.exit("Unexpected number of columns in the DOSCAR total DOS block.")

    dos = np.array(dos)
    energy = np.append(dos[:, 0], dos[-1, 0])
    energy = np.insert(energy, 0, energy[0])
    dos_up = close_polygon(dos[:, 1])
    dos_down = close_polygon(dos[:, 2]) if is_magnetic else None
    return energy, dos_up, dos_down, e_fermi, is_magnetic


def read_eigenvalues(eigenfile, is_magnetic):
    """Return (bands_up, bands_down, nbands). bands_down is None if not spin-polarized."""
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
    bands_up, bands_down = [], []
    for j in range(nbands):
        rows = [lines[i].split() for i in range(j, nkpt * block, block)]
        bands_up.append([float(r[1]) for r in rows])
        if is_magnetic:
            bands_down.append([float(r[2]) for r in rows])
    return bands_up, (bands_down if is_magnetic else None), nbands


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
    names.append(contents[-1].split()[-1].strip('!#'))
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
    plt.rc('font', serif=['Times New Roman', 'DejaVu Serif'])
    plt.rcParams['axes.linewidth'] = 1.20
    for axis in ('x', 'y'):
        plt.rcParams[axis + 'tick.major.size'] = 6
        plt.rcParams[axis + 'tick.minor.size'] = 3
        plt.rcParams[axis + 'tick.major.width'] = 1.20
        plt.rcParams[axis + 'tick.minor.width'] = 1.20
        plt.rcParams[axis + 'tick.direction'] = 'in'


def title_axes(fig, title):
    """Invisible axes spanning the figure, used to place a centered title."""
    ax = fig.add_subplot(111)
    for side in ('top', 'bottom', 'left', 'right'):
        ax.spines[side].set_color('none')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    return ax


def draw_bands(ax, bands, kdist, kdist_hsym, kname, fermi, ymin, ymax, color=None):
    kmax = max(kdist)
    for band in bands:
        kwargs = {'lw': 1}
        if color:
            kwargs['c'] = color
        ax.plot(kdist, np.array(band) - fermi, **kwargs)
    ax.axhline(0, color='black', lw=0.5, ls='--')
    ax.axis([0, kmax, ymin, ymax])
    for x in kdist_hsym:
        ax.axvline(x, color='k', lw=0.5, ls='--')
    ax.set_xticks(kdist_hsym)
    ax.set_xticklabels(kname)


def main():
    args = parse_args()
    check_inputs()

    energy, dos_up, dos_down, fermi_doscar, is_magnetic = read_total_dos('DOSCAR')
    bands_up, bands_down, nbands = read_eigenvalues('EIGENVAL', is_magnetic)
    kname, points_per_segment = read_kpath_labels('KPOINTS')
    fermi_outcar, kpoints = read_outcar('OUTCAR')
    kdist, kdist_hsym = k_distances(kpoints, points_per_segment)

    source = args.fermi
    if source == 'auto':
        source = 'doscar' if is_magnetic else 'outcar'
    fermi = fermi_doscar if source == 'doscar' else fermi_outcar
    print("Fermi level (%s): %.4f eV" % (source, fermi))

    ymin = args.emin if args.emin is not None else (-4.0 if is_magnetic else -3.0)
    ymax = args.emax if args.emax is not None else (4.0 if is_magnetic else 3.0)
    dosmax = args.dosmax

    set_style()
    fig = plt.figure()
    rel = energy - fermi

    if is_magnetic:
        gs = gridspec.GridSpec(1, 3, width_ratios=[3, 3, 1.2])
        ax = title_axes(fig, args.title)
        ax0 = fig.add_subplot(gs[0])   # spin down bands
        ax1 = fig.add_subplot(gs[1])   # spin up bands
        ax2 = fig.add_subplot(gs[2])   # DOS

        draw_bands(ax1, bands_up, kdist, kdist_hsym, kname, fermi, ymin, ymax, 'blue')
        draw_bands(ax0, bands_down, kdist, kdist_hsym, kname, fermi, ymin, ymax, 'red')
        ax0.set_ylabel(r"$E-E_{\mathrm{F}}$ [eV]", fontsize=12)

        ax2.plot(dos_up, rel, 'b-', lw=0.5)
        ax2.plot(-dos_down, rel, 'r-', lw=0.5)
        ax2.fill(dos_up, rel, lw=0, facecolor='b', alpha=0.8)
        ax2.fill(-dos_down, rel, lw=0, facecolor='r', alpha=0.8)
        ax2.axhline(0, color='black', lw=0.5, ls='--')
        ax2.axis([-dosmax, dosmax, ymin, ymax])
        ax2.set_xticks(np.array([-0.6, 0.0, 0.6]) * dosmax)
        ax2.axvline(0, color='k')
        wspace = 0.05
        hide_y = (ax1, ax2)
    else:
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1])
        ax = title_axes(fig, args.title)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])

        draw_bands(ax1, bands_up, kdist, kdist_hsym, kname, fermi, ymin, ymax)
        ax1.set_ylabel(r"$E-E_{\mathrm{F}}$ [eV]")

        ax2.plot(dos_up, rel, 'b-', lw=0.5)
        ax2.axhline(0, color='black', lw=0.5, ls='--')
        ax2.axis([0.0, dosmax, ymin, ymax])
        ax2.set_xticks(np.arange(0.1, 1.0, 0.2) * dosmax)
        ax2.fill(dos_up, rel, lw=0, facecolor='b', alpha=0.8)
        ax2.axvline(0, color='k')
        wspace = 0.03
        hide_y = (ax2,)

    for a in hide_y:
        a.tick_params(labelleft=False)

    plt.subplots_adjust(left=0.135, bottom=0.06, top=0.921, right=0.935,
                        hspace=0.0, wspace=wspace)
    plt.savefig(args.outname + '.pdf', bbox_inches='tight', pad_inches=0.01)
    plt.close(fig)


if __name__ == '__main__':
    main()

