"""
Run crystals in energy_vs_volume
"""

# Find experimental... OR better, DFT data on equations of state for cubic cell structures
# MgO, NaCl, Si, Diamond, Ge (expect to break), TiO2 (potentially), Copper
# Also want experimental and DFT lattice constants

# Fit the Birch-Murnaghan equation to each pair of data

from energy_vs_volume import ev_functions

ev_functions.delta_coefficients([1,2,3,4])
