"""   Compare the energy per atom for system described by its primitive unit cell
      and its conventional unit cell:
      The energy per atom should be the same, independent of the choice of unit cell.

      FCC systems are good candidates as it's simple to move from the 2-atom primitive
      unit cell to the 8-atom conventional unit cell.
"""

from collections import OrderedDict

from crystal_system import cubic
from src.pymatgen_wrappers import cif_parser_wrapper
from src.xtb_potential import xtb_potential_str, PotentialType
from src import qcore_input_strings as qcore_input
from src.utils import Set, SetAssert


# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED


def primitive_silicon(named_result='primitive_cell_si') -> str:
    """
    Primitive silicon

    Note, this will likely give slightly difference results to using
    cubic.silicon, which (probably) has a slightly different lattice
    parameter

    :param named_result: named_result
    :return: input string
    """

    fname = '../' + cubic.conventional_fcc_cifs['silicon'].file
    crystal = cif_parser_wrapper(fname, is_primitive_cell=True, fractional=True,
                                 bravais='fcc')
    assert crystal['n_atoms'] == 2, "Expect 2-atom basis for primitive si"

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return qcore_input.xtb_input_string_cleaner(crystal=crystal, options=options,
                                                assertions=assertions[potential_type],
                                                named_result=named_result)



def conventional_silicon(named_result='conventional_cell_si'):

    fname = '../' + cubic.conventional_fcc_cifs['silicon'].file
    crystal = cif_parser_wrapper(fname, is_primitive_cell=False, fractional=True,
                                 bravais='cubic')
    assert crystal['n_atoms'] == 8, "Expect 8-atom basis for primitive si"

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    return qcore_input.xtb_input_string_cleaner(crystal=crystal, options=options,
                                                named_result=named_result)


def primitive_copper():

    named_result = 'primitive_copper'
    fname = '../' + cubic.fcc_cifs['copper'].file
    crystal = cif_parser_wrapper(fname,  is_primitive_cell=True, fractional=True)
    settings = OrderedDict([
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
    settings = OrderedDict([
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
    settings = OrderedDict([
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
    settings = OrderedDict([
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