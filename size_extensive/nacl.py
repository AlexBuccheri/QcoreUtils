""" Sodium chloride conventional cell for 3 supercells"""

from collections import OrderedDict

from src.xtb_potential import xtb_potential_str, PotentialType
from src.utils import Set, SetAssert
from crystal_system import cubic
from src.pymatgen_wrappers import cif_parser_wrapper
from src.qcore_input_strings import xtb_input_string_cleaner


# Use correct lattice vectors
def get_supercell_bravais(supercell_coefficients: list):
    """
    Input assumes a simple cubic cell, however if one does not expand
    isotropically, then the resulting super cell may become
    simple tetragonal (i.e. [n, n, m]) or
    simple orthorhombic (i.e. [n, m, l])
    """
    assert len(supercell_coefficients) == 3, "len(supercell_coefficients) != 3"
    n_unique = len(set(supercell_coefficients))
    if n_unique == 1:
        return 'cubic'
    elif n_unique == 2:
        return 'tetragonal'
    elif n_unique == 3:
        return 'orthorhombic'


# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED


def convention_sodium_chloride(supercell_coefficients, named_result):

    fname = '../' + cubic.conventional_fcc_cifs['sodium_chloride'].file
    crystal = cif_parser_wrapper(fname,
                                 is_primitive_cell=False,
                                 fractional=True,
                                 bravais=get_supercell_bravais(supercell_coefficients),
                                 supercell_coefficients=supercell_coefficients)

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([4, 4, 4])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin')),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return xtb_input_string_cleaner(crystal=crystal, options=options,
                                    assertions=assertions[potential_type],
                                    named_result=named_result)


conventional_cell = convention_sodium_chloride([1, 1, 1], 'conv_cell')

cubic_super_cell = convention_sodium_chloride([2, 2, 2], 'supercell')

print(conventional_cell)

