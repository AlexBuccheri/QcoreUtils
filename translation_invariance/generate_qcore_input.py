"""
Generate translational invariance input file.

See translation_invariance.crystals for each crystal setting
and details on the app test.
"""
from translation_invariance.crystals import *

def produce_translation_test(crystal_inputs):
    total_input = ''
    for input in crystal_inputs:
        total_input += input() + '\n'
    return total_input


test_header = \
"""! Test translational invariance of the total energy for periodic GFN1-xTB.
! One translational invariance test per Bravais lattice
! For each material, do two calculations and assert that the energies are the same

! The shifts in positions are arbitrary. If a shift results in
! any fractional position exceeding 1, 'wrap_atoms = true' should be used

! TODO(Alex) Add simple, face-centred & base-centred orthorhombic tests
! TODO(Alex) Add simple, and base-centred C monoclinic tests
! TODO(Alex) Add triclinic lattice test
"""

# TODO(Alex) Add remaining lattices
all_tests = [alpha_N2,               # Simple cubic
             sio2_zeolite,           # BCC
             silicon,                # FCC
             copper,                 # FCC metal
             tio2_rutile,            # Simple tetragonal
             tio2_anatase,           # Base-centred tetragonal
             boron_nitride_hex,      # Hexagonal
             molybdenum_disulfide,   # Rhombohedral
             copper_oxyfluoride      # Body-centred orthorhombic
             ]

print(test_header)
print(produce_translation_test(all_tests))
#print(produce_translation_test([alpha_N2]))

