"""   Compare the energy per atom for system described by its primitive unit cell
      and its conventional unit cell:
      The energy per atom should be the same, independent of the choice of unit cell.
"""

import collections

from crystal_system import cubic
from src.pymatgen_wrappers import cif_parser_wrapper
from src import qcore_input_strings as qcore_input
from src.utils import Set


def primitive_copper():

    named_result = 'primitive_copper'
    fname = '../' + cubic.fcc_cifs['copper'].file
    crystal = cif_parser_wrapper(fname,  is_primitive_cell=True, fractional=True)
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([6, 6, 6])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    input_string = qcore_input.xtb_input_string(crystal, settings, named_result=named_result)
    return input_string

def conventional_copper():

    named_result = 'conventional_copper'
    fname = '../' + cubic.conventional_fcc_cifs['copper'].file
    crystal = cif_parser_wrapper(fname, is_primitive_cell=False, fractional=True, bravais='cubic')
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([4, 4, 4])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    input_string = qcore_input.xtb_input_string(crystal, settings, named_result=named_result)
    return input_string


def primitive_sodium_chloride():

    named_result = 'primitive_nacl'
    fname = '../' + cubic.fcc_cifs['sodium_chloride'].file
    crystal = cif_parser_wrapper(fname,  is_primitive_cell=True, fractional=True)
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([4, 4, 4])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    input_string = qcore_input.xtb_input_string(crystal, settings, named_result=named_result)
    return input_string

def convention_sodium_chloride():

    named_result = 'conventional_nacl'
    fname = '../' + cubic.conventional_fcc_cifs['sodium_chloride'].file
    crystal = cif_parser_wrapper(fname, is_primitive_cell=False, fractional=True, bravais='cubic')
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([4, 4, 4])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    input_string = qcore_input.xtb_input_string(crystal, settings, named_result=named_result)
    return input_string




input_copper = primitive_copper() + '\n\n' + conventional_copper()
input_nacl = primitive_sodium_chloride() + '\n\n' + convention_sodium_chloride()

print(input_nacl)


