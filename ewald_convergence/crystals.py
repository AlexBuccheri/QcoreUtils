"""
Inputs for checking Ewald convergence
"""

import collections

from crystal_system import cubic
from src.pymatgen_wrappers import cif_parser_wrapper
from src import qcore_input_strings as qcore_input
from src.utils import Set

def convention_magnesium_oxide(named_result, ewald):
    fname = '../' + cubic.conventional_fcc_cifs['magnesium_oxide'].file
    crystal = cif_parser_wrapper(fname, is_primitive_cell=False, fractional=True, bravais='cubic')
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(ewald['real'], 'bohr')),
        ('ewald_reciprocal_cutoff', Set(ewald['reciprocal'])),
        ('ewald_alpha',             Set(ewald['alpha'])),
        ('monkhorst_pack',          Set([8, 8, 8])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    input_string = qcore_input.xtb_input_string(crystal, settings, named_result=named_result)
    return input_string
