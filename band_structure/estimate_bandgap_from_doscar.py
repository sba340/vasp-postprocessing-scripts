#!/usr/bin/env python3
"""Estimate band gaps from a VASP DOSCAR file.

Run this script in a calculation directory that contains DOSCAR and OSZICAR.
It finds the energy grid point closest to the Fermi level, then walks along
the grid while the density of states is zero. It prints one line with the
directory name, a value read from OSZICAR, two gap estimates in eV, the DOS
columns at the Fermi level, and a flag (I = insulator, M = metal).

Usage:
    python estimate_bandgap_from_doscar.py
"""
import os

dirname = os.path.basename(os.getcwd())

with open('OSZICAR') as f:
    osfile = f.readlines()

dos = []
with open('DOSCAR') as fdos:
    for i, line in enumerate(fdos, start=1):
        if i == 6:
            _emax, _emin, nedos, eferm, _weight = line.split()
            nedos = float(nedos)
            eferm = float(eferm)
        if 6 < i < nedos + 6:
            dos.append([float(x) for x in line.split()])

# Grid point closest to the Fermi level
i_eferm = 0
del_e = abs(dos[0][0] - eferm)
for i in range(len(dos)):
    if abs(dos[i][0] - eferm) < del_e:
        del_e = abs(dos[i][0] - eferm)
        i_eferm = i

# Walk outward while the DOS is zero
n = i_eferm
while dos[n + 1][1] == 0:
    n += 1

m = i_eferm
while dos[m + 1][2] == 0:
    m += 1

k = i_eferm
while dos[k - 1][2] == 0:
    k -= 1

bandgap1 = dos[n][0] - dos[i_eferm][0]
bandgap2 = dos[m][0] - dos[k][0]
state = 'I' if bandgap1 > 0.0 else 'M'

print('%12s %12s %10s %8.4f %8.4f %8.4f %8.4f %6s' % (
    dirname, osfile[-2].split()[3], 'gap:  ', bandgap1, bandgap2,
    dos[i_eferm][1], dos[i_eferm][2], state))

