
#!/usr/bin/env python3
"""Plot the total density of states from a VASP DOSCAR file.

Run this script in a directory that contains DOSCAR. The script detects
whether the calculation is spin-polarized from the number of columns in the
total DOS block (3 columns: non-spin-polarized; 5 columns: spin-polarized).
Energies are plotted relative to the Fermi level. For spin-polarized runs the
spin-down DOS is drawn as negative values.

Output:
    vasp_dos_nonmagnetic.pdf or vasp_dos_magnetic.pdf (or <outname>.pdf)

Usage:
    python plot_total_dos_from_doscar.py -emin -1 -emax 1
"""
import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument('-emin', dest='emin', type=float, default=-1.0,
                        help='Lower limit of the energy axis relative to E_F '
                             '(default: -1)')
    parser.add_argument('-emax', dest='emax', type=float, default=1.0,
                        help='Upper limit of the energy axis relative to E_F '
                             '(default: 1)')
    parser.add_argument('-outname', dest='outname', type=str, default=None,
                        help='Output file name without extension')
    return parser.parse_args()


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


def set_style():
    plt.rc('font', family='serif')
    plt.rc('font', serif=['Times New Roman', 'DejaVu Serif'], size=11)
    plt.rcParams['axes.linewidth'] = 0.7
    for axis in ('x', 'y'):
        plt.rcParams[axis + 'tick.major.size'] = 4
        plt.rcParams[axis + 'tick.minor.size'] = 2
        plt.rcParams[axis + 'tick.major.width'] = 0.7
        plt.rcParams[axis + 'tick.minor.width'] = 0.7


def main():
    args = parse_args()
    energy, dos_up, dos_down, e_fermi, is_magnetic = read_total_dos('DOSCAR')

    rel = energy - e_fermi
    window = (rel > args.emin) & (rel < args.emax)

    set_style()
    fig = plt.figure(figsize=(4, 2.5))
    ax = fig.add_subplot(111)

    ax.plot(rel, dos_up, lw=0.5, color="#F7A416")
    ax.fill(rel, dos_up, lw=0, facecolor='#F7CB15', alpha=0.8)

    if is_magnetic:
        ax.plot(rel, -dos_down, lw=0.5, color='#243010')
        ax.fill(rel, -dos_down, lw=0, facecolor='#81C14B', alpha=0.8)
        ax.set_title('Magnetic calculation')
        ax.set_ylim([-np.amax(dos_down[window]), np.amax(dos_up[window])])
        default_name = 'vasp_dos_magnetic'
    else:
        ax.set_title('Non-magnetic calculation')
        ax.set_ylim([0, np.amax(dos_up[window])])
        default_name = 'vasp_dos_nonmagnetic'

    ax.set_xlabel(r"$E-E_{\mathrm{F}}$ [eV]")
    ax.set_ylabel(r"DOS [states eV$^{-1}$ cell$^{-1}$]")
    ax.axvline(0, color='black', lw=0.5, ls='--')
    ax.set_xlim([args.emin, args.emax])

    plt.savefig((args.outname or default_name) + '.pdf',
                bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)


if __name__ == '__main__':
    main()
