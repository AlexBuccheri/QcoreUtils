"""
Inputs  for asserting symmetry_reduction of the Monkhorst Pack grid
Currently restricted to cubic lattices
"""

from collections import OrderedDict

from .app_test_string import xtb_symmetry_test_string

from crystal_system import cubic, hexagonal
from src.xtb_potential import xtb_potential_str, PotentialType
from src.pymatgen_wrappers import cif_parser_wrapper
from src.utils import Set, SetAssert


# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED


def alpha_N2(named_result='alpha_N2') -> str:
    """
    Simple cubic
    :param named_result: named result string
    :return: Input for testing symmetry reduction of MP grid
    """

    fname = '../' + cubic.cubic_cifs['alpha-N2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'cubic', "crystal not simple cubic"

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin')),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return xtb_symmetry_test_string(crystal=crystal,
                                    options=options,
                                    assertions=assertions[PotentialType.TRUNCATED],
                                    named_result=named_result)

def potassium(named_result='potassium') -> str:
    """
     BCC
    :param named_result: named result string
    :return: Input for testing symmetry reduction of MP grid
    """

    fname = '../' + cubic.bcc_cifs['potassium'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not bcc"

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin')),
        ('solver',                  Set('SCC')),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return xtb_symmetry_test_string(crystal=crystal,
                                    options=options,
                                    assertions=assertions[PotentialType.TRUNCATED],
                                    named_result=named_result)


def nickel_triphosphide(named_result='nickel_triphosphide') -> str:
    """
     BCC
     Won't converge with these settings
    :param named_result: named result string
    :return: Input for testing symmetry reduction of MP grid
    """

    fname = '../' + cubic.bcc_cifs[named_result].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not bcc"

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin')),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return xtb_symmetry_test_string(crystal=crystal,
                                    options=options,
                                    assertions=assertions[PotentialType.TRUNCATED],
                                    named_result=named_result)

def sio2_zeolite(named_result='SiO2') -> str:
    """
    BCC structure
    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    arbitrary_shift = 0.1

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(40, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('potential_type',          Set(xtb_potential_str(potential_type))),
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(14)),
                                              ("energy", SetAssert(-69.074163659, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    comments = '! Fractional shift of 0.1'

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result, comments=comments)


def silicon(named_result='silicon') -> str:
    """
     FCC
    :param named_result: named result string
    :return: Input for testing symmetry reduction of MP grid
    """

    crystal = cubic.silicon()
    assert crystal['bravais'] == 'fcc', "crystal not fcc"

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin')),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return xtb_symmetry_test_string(crystal=crystal,
                                    options=options,
                                    assertions=assertions[PotentialType.TRUNCATED],
                                    named_result=named_result)


# Hexagonal
def boron_nitride_hex(named_result='boron_nitride') -> str:
    """
     Hexagonal
    :param named_result: named result string
    :return: Input for testing symmetry reduction of MP grid
    """

    crystal = hexagonal.boron_nitride()
    assert crystal['bravais'] == 'hexagonal', "crystal not hexagonal"

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin')),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return xtb_symmetry_test_string(crystal=crystal,
                                    options=options,
                                    assertions=assertions[PotentialType.TRUNCATED],
                                    named_result=named_result)



