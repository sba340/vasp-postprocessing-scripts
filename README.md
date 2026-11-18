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
| quick_band_plot_and_gap_from_eigenval.py | band_structure | Quick band plot by k-point index with band edge and gap summary from EIGENVAL and OUTCAR (non-spin-polarized) |
| plot_total_dos_from_doscar.py | band_structure | Plots the total density of states from DOSCAR, spin-polarized or not |
| plot_band_structure_with_dos.py | band_structure | Plots a band structure next to the total DOS from EIGENVAL, KPOINTS, OUTCAR and DOSCAR, spin-polarized or not |

Input files are not included. Each script reads standard VASP output from the
current working directory.

