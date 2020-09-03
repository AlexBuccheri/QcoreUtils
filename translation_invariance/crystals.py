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
# TODO(Alex) This has a lot of atoms with positions at 1 => MUST ADD WRAPPING TO ENTOS
# TODO(Alex) Need to set assertion values by running first
"""


from collections import OrderedDict
import enum

from src.utils import Set, SetAssert
from crystal_system import cubic, tetragonal, hexagonal, trigonal
from src.pymatgen_wrappers import cif_parser_wrapper
from translation_invariance.string_generator import xtb_translational_invariance_string as translation_string


class PotentialType(enum.Enum):
    TRUNCATED = enum.auto
    FULL = enum.auto


# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED


def set_xtb_potential(potential_type: enum.Enum) -> dict:
    """
    Set xTB potential
    :param potential_type:
    :return:
    """
    if potential_type == PotentialType.TRUNCATED:
        return OrderedDict([
        ('potential_type', Set('truncated')),
        ('smoothing_range', Set(1, 'bohr'))
    ])
    elif potential_type == PotentialType.Full:
        return OrderedDict([('potential_type', Set('full'))])
    else:
        raise Exception("Invalid choice for potential_type: ", potential_type)


def alpha_N2(named_result='alpha_N2') -> str:
    """
    Simple cubic
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    fname = '../' + cubic.cubic_cifs['alpha-N2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'cubic', "crystal not simple cubic"

    # Can't do a larger shift as there's no wrapping of atomic positions in entos
    arbitrary_shift = 0.05

    options = OrderedDict([
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

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0, 0)), ("energy", SetAssert(0, 0))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)), ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, sub_commands, options, assertions[potential_type],
                              arbitrary_shift, named_result)


def malh_perovskite(named_result='malh_perovskite'):
    """
    Second simple cubic, with several atoms in the central cell.
    #TODO Relaxed structure looks triclinic from lattice parameters

    Methylammonium lead halides (MALHs)
    chemical formula of CH3NH3PbX3, where X = I, Br or Cl.
    CH3NH3PbI3

    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + cubic.cubic_cifs['CH3NH3PbI3'].file
    crystal = cif_parser_wrapper(fname)
    # Relaxed structure looks triclinic from lattice parameters
    assert crystal['bravais'] == 'cubic', "crystal not simple cubic"

    #TODO(Alex) Set an arb shift once wrapping of positions can be done
    arbitrary_shift = 0.0

    options = OrderedDict([
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

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0, 0)), ("energy", SetAssert(0, 0))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)), ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, sub_commands, options, assertions[potential_type],
                              arbitrary_shift, named_result)



def sio2_zeolite(named_result='SiO2') -> str:
    """
    BCC structure
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    # TODO(Alex) This has a lot of atoms with positions at 1 => MUST ADD WRAPPING
    arbitrary_shift = 0.0

    options = OrderedDict([
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

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0, 0)), ("energy", SetAssert(0, 0))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)), ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, sub_commands, options, assertions[potential_type],
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
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0, 0)), ("energy", SetAssert(0, 0))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)), ("energy", SetAssert(0, 0))])
    }


    return translation_string(crystal, sub_commands, options, assertions[potential_type],
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
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter",  SetAssert(3)),
                                              ("energy",  SetAssert(-3.827222, 1.e-6)),
                                              ("entropy", SetAssert(-0.0001098, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(3)),
                                              ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, sub_commands, options, assertions[potential_type],
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
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(35)),
                                                        ("energy", SetAssert(-150.135037, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, sub_commands, options, assertions[potential_type],
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
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(11)),
                                                        ("energy", SetAssert(-19.142514, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, sub_commands, options, assertions[potential_type],
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
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, sub_commands, options, assertions[potential_type],
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
        ('ewald_reciprocal_cutoff', Set(2)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    sub_commands = OrderedDict([('xtb_potential', set_xtb_potential(potential_type))])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, sub_commands, options, assertions[potential_type],
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

