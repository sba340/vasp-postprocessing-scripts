# VASP Post-processing Scripts

Python scripts for analyzing and plotting output from VASP density functional
theory calculations: band gaps, band structures, density of states, bonding
analysis, and phonons.

## Requirements

Python 3.8 or later. Individual scripts list their own dependencies.

## Scripts

| Script | Folder | Description |
|---|---|---|
| estimate_bandgap_from_doscar.py | band_structure | Estimates band gaps from DOSCAR |
| plot_band_structure_from_eigenval.py | band_structure | Plots a band structure from EIGENVAL, KPOINTS and OUTCAR (non-spin-polarized) |

Input files are not included. Each script reads standard VASP output from the
current working directory.
