"""
Generate translational invariance input file.

See translation_invariance.crystals for each crystal setting
and details on the app test.
"""
from translation_invariance.crystals import *

# TODO Add remaining 3 lattices: monoclinic and triclinic
all_tests = [alpha_N2,               # Simple cubic
             sio2_zeolite,           # BCC
             silicon,                # FCC
             copper,                 # FCC metal (superfluous for this test)
             tio2_rutile,            # Simple tetragonal
             tio2_anatase,           # Base-centred tetragonal
             boron_nitride_hex,      # Hexagonal
             molybdenum_disulfide,   # Rhombohedral
             gold_cadmium,           # Simple orthorhombic
             copper_oxyfluoride,     # Body-centred orthorhombic
             aluminium_titanate,     # Base-centred C orthorhombic
             titanium_disilicide     # Face-centred orthorhombic.
             ]


def produce_translation_test(crystal_inputs):
    total_input = ''
    for input in crystal_inputs:
        total_input += input() + '\n'
    return total_input


print(produce_translation_test([titanium_disilicide]))



