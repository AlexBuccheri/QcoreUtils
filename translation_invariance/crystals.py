"""
Test translational invariance for periodic GFN1 xTB in ENTOS Qcore.
One translational invariance test per Bravais lattice

All inputs are specified in case defaults change in the future:
Better to be explicit

Options are the same in most cases - we're not interested in converged calculations.
They simply need to be fast, hence default (H,S, rep) cut-offs, small Ewald cut-offs
and small MP grid density. Some will use SCC or symmetry for faster calculations.

Potential type determines how the real-space sum in the 2nd-order xTB potential is
handled. One can switch between potential_type with the variable at the top of the script.

arbitrary_shifts are given in fractional units [0:1]
"""

from collections import OrderedDict
import enum

from src.utils import Set, SetAssert
from crystal_system import cubic, tetragonal, hexagonal, trigonal
from src.pymatgen_wrappers import cif_parser_wrapper
from translation_invariance.string_generator import xtb_translational_invariance_string as translation_string


def xtb_potential_str(potential_type: enum.Enum) -> str:
    """
    Set xTB potential
    :param potential_type:
    :return:
    """
    if potential_type == PotentialType.TRUNCATED:
        return 'truncated'

    elif potential_type == PotentialType.Full:
        return 'full'
    else:
        raise Exception("Invalid choice for potential_type: ", potential_type)


class PotentialType(enum.Enum):
    TRUNCATED = 1
    FULL = 2


# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED


def alpha_N2(named_result='alpha_N2') -> str:
    """
    Simple cubic
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    fname = '../' + cubic.cubic_cifs['alpha-N2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'cubic', "crystal not simple cubic"
    arbitrary_shift = 0.05

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
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(2)),
                                             ("energy", SetAssert(-24.734788141, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, options, assertions[PotentialType.TRUNCATED],
                              arbitrary_shift, named_result)


def malh_perovskite(named_result='malh_perovskite'):
    """
    Second simple cubic, with several atoms in the central cell.
    #TODO Relaxed structure looks triclinic from lattice parameters
    # Need to modify alpha to 89.9 but it does converge

    Methylammonium lead halides (MALHs)
    chemical formula of CH3NH3PbX3, where X = I, Br or Cl.
    CH3NH3PbI3

    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + cubic.cubic_cifs['CH3NH3PbI3'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'monoclinic', "crystal not simple monoclinic"

    arbitrary_shift = 0.65

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
        ('potential_type',          Set(xtb_potential_str(potential_type))),
        ('wraps_atoms',             Set(True))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)


def sio2_zeolite(named_result='SiO2') -> str:
    """
    BCC structure
    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    #TODO(Alex) Add ('wraps_atoms', Set(True)) to structure command

    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    arbitrary_shift = 0.51

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
        ('potential_type', Set(xtb_potential_str(potential_type))),
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(16)),
                                              ("energy", SetAssert(-68.109016142, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)



def silicon(named_result='silicon') -> str:
    """
    FCC structure
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    crystal = cubic.silicon()
    arbitrary_shift = 0.43

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
        ('potential_type', Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(8)),
                                              ("energy", SetAssert(-3.669583993, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)


def copper(named_result='copper') -> str:
    """
    FCC metal
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    fname = '../' + cubic.fcc_cifs['copper'].file
    crystal = cif_parser_wrapper(fname, fractional=True)
    arbitrary_shift = 0.53

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
        ('potential_type', Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter",  SetAssert(3)),
                                              ("energy",  SetAssert(-3.080165805, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(3)),
                                              ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)


def tio2_rutile(named_result='rutile') -> str:
    """
    Simple tetragonal
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    crystal = tetragonal.tio2_rutile()
    arbitrary_shift = 0.134

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin')),
        ('solver',                  Set('SCC')),
        ('potential_type', Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(14)),
                                                        ("energy", SetAssert(-22.337682356, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)


def tio2_anatase(named_result='anatase') -> str:
    """
    Base-centred tetragonal
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    crystal = tetragonal.tio2_anatase()
    arbitrary_shift = 0.04

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin')),
        ('potential_type', Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(11)),
                                                        ("energy", SetAssert(-19.256467094, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)


def boron_nitride_hex(named_result='bn_hex') -> str:
    """
    Hexagonal.
    Check D3 is added
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    crystal = hexagonal.boron_nitride()
    arbitrary_shift = 0.13

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin')),
        ('potential_type', Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(7)),
                                                        ("energy", SetAssert(-9.426842352, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)


def molybdenum_disulfide(named_result='MoS2') -> str:
    """
    Rhombohedral.
    Check D3 is added
    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    crystal = trigonal.MoS2()
    arbitrary_shift = 0.53

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(10, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin')),
        ('potential_type', Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(36)),
                                                        ("energy", SetAssert(-76.45562556, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)



# #
# # TODO Add simple orthorhombic
#
# # TODO Find: base-centred C orthorhombic
#
# # TODO Find: body-centred orthorhombic
#
# # TODO Find: face-centred orthorhombic
#
# # TODO Find: simple monoclinic
#
# # TODO Find: base-centred C monoclinic
#
#
# def produce_translation_test(crystal_inputs):
#     total_input = ''
#     for input in crystal_inputs:
#         total_input += input() + '\n'
#     return total_input
#
#
# all_tests = [silicon, copper, tio2_rutile, tio2_anatase, boron_nitride_hex, molybdenum_disulfide_rhom]
# print(produce_translation_test([silicon]))

