""" qCore inputs for crystals, testing translational invariance """

import collections

from utils import Set
from crystal_system import cubic, tetragonal, hexagonal, trigonal, orthorhombic, monoclinic, triclinic
from pymatgen_wrappers import cif_parser_wrapper
from translation_invariance.string_generator import xtb_translational_invariance_string as translation_string


# TODO Add: simple cubic

# TODO Add: bcc

# fcc
def silicon() -> str:

    crystal = cubic.silicon()
    named_result = 'silicon'
    arbitrary_shift = 0.43

    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    assertions = collections.OrderedDict([
        ("n_iter", 8),
        ("energy", -3.669686)
    ])

    return translation_string(crystal, settings, assertions, arbitrary_shift, named_result)


# fcc, metal
def copper() -> str:

    fname = '../' + cubic.fcc_cifs['copper'].file
    crystal = cif_parser_wrapper(fname, fractional=True)

    named_result = 'copper'
    arbitrary_shift = 0.53

    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    assertions = collections.OrderedDict([
        ("n_iter",   3),
        ("energy",  -3.827222),
        ("entropy", -0.0001098)
        ])

    return translation_string(crystal, settings, assertions, arbitrary_shift, named_result)


# simple tetragonal
def tio2_rutile() -> str:

    crystal = tetragonal.tio2_rutile()
    named_result = 'rutile'
    arbitrary_shift = 0.134

    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    assertions = collections.OrderedDict([
        ("n_iter", 35),
        ("energy", -150.135037)
    ])

    return translation_string(crystal, settings, assertions, arbitrary_shift, named_result)


# body-centred tetragonal
def tio2_anatase() -> str:

    crystal = tetragonal.tio2_anatase()
    named_result = 'anatase'
    arbitrary_shift = 0.04

    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    assertions = collections.OrderedDict([
        ("n_iter", 11),
        ("energy", -19.142514)
    ])

    return translation_string(crystal, settings, assertions, arbitrary_shift, named_result)

# hexagonal
def boron_nitride_hex() -> str:

    crystal = hexagonal.boron_nitride()
    named_result = 'bn_hex'
    comments = "! Periodic D3 should be included for an accurate total energy"
    arbitrary_shift = 0.13

    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    assertions = collections.OrderedDict([
        ("n_iter", 6),
        ("energy", -9.431841)
    ])

    return translation_string(crystal, settings, assertions, arbitrary_shift, named_result, comments)

# rhombohedral
def molybdenum_disulfide_rhom() -> str:

    crystal = trigonal.MoS2()
    named_result = 'MoS2_rhom'
    comments = "! Periodic D3 should be included for an accurate total energy"
    arbitrary_shift = 0.53

    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    assertions = collections.OrderedDict([
        ("n_iter", 24),
        ("energy", -9.505776937)
    ])

    return translation_string(crystal, settings, assertions, arbitrary_shift, named_result, comments)

# simple orthorhombic
# TODO Try one of Rui's

# TODO Find: base-centred C orthorhombic

# TODO Find: body-centred orthorhombic

# TODO Find: face-centred orthorhombic

# TODO Find: simple monoclinic

# TODO Find: base-centred C monoclinic

# TODO Find: triclinic



def produce_translation_test(crystal_inputs):
    total_input = ''
    for input in crystal_inputs:
        total_input += input() + '\n'
    return total_input


all_tests = [silicon, copper, tio2_rutile, tio2_anatase, boron_nitride_hex, molybdenum_disulfide_rhom]
print(produce_translation_test([silicon]))

