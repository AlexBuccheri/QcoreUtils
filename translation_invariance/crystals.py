"""
Test translational invariance for periodic GFN1 xTB in ENTOS Qcore.
One translational invariance test per Bravais lattice

All inputs are specified in case defaults change in the future:
Better to be explicit

Options are the same in most cases - we're not interested in converged calculations.
They simply need to be fast, hence default (H,S, rep) cut-offs, small Ewald cut-offs
and small MP grid density. Some will use SCC or symmetry (where possible) for faster
calculations.

Potential type determines how the real-space sum in the 2nd-order xTB potential is
handled. One can switch between potential_type with the variable at the top of the script.

arbitrary_shifts are given in fractional units and can be anything in principle.
If A shifted position lies outside of [0:1], then wrap_atoms = true is automatically
added to the structure command
"""

from collections import OrderedDict

from src.xtb_potential import xtb_potential_str, PotentialType
from src.utils import Set, SetAssert
from crystal_system import cubic, tetragonal, hexagonal, trigonal, orthorhombic
from src.pymatgen_wrappers import cif_parser_wrapper
from translation_invariance.string_generator \
    import xtb_translational_invariance_string as translation_string, \
    comments_list


# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED


def alpha_N2(named_result='alpha_N2') -> str:
    """
    Simple cubic
    TODO(Alex) If the shift is large for alpha N2, the total energy is massively different
    (~0.5 bohr). Also, the enthropy is zero for the translated system:
    Looks like an issue with occupations. Requires debugging at Gamma

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
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(2)),
                                             ("energy", SetAssert(-24.742649248, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result, comments=comments_list(arbitrary_shift))


def sio2_zeolite(named_result='SiO2') -> str:
    """
    BCC structure
    Notes:
    If one exceeds 0.24, one has to greatly increase the cut-offs to
    get consistent total energies.
    See range_of_translations.py for more exploration with this system

    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    arbitrary_shift = 0.24

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
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(17)),
                                              ("energy", SetAssert(-69.056517465, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result, comments=comments_list(arbitrary_shift))


def silicon(named_result='silicon') -> str:
    """
    FCC structure
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    crystal = cubic.silicon()
    arbitrary_shift = 0.5

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
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(7)),
                                              ("energy", SetAssert(-3.669583993, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                              ("energy", SetAssert(0, 0))])
    }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result, comments=comments_list(arbitrary_shift))


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
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter",  SetAssert(3)),
                                              ("energy",  SetAssert(-3.080495065, 1.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(3)),
                                              ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result, comments=comments_list(arbitrary_shift))


def tio2_rutile(named_result='rutile') -> str:
    """
    Simple tetragonal
    Notes:
    Similar issue to SiO2. If one goes to a larger translation where
    wrapping is required (0.24), then one has to turn up the cut-offs
    to obtain agreement to 7 d.p. in the total_energy.
    See the gamma_point input in range_of_translations.py

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
        ('solver',                  Set('SCC')),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(14)),
                                                        ("energy", SetAssert(-22.336229856, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


def tio2_anatase(named_result='anatase') -> str:
    """
    Base-centred tetragonal
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    crystal = tetragonal.tio2_anatase()
    arbitrary_shift = 0.34

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(40, 'bohr')),
        ('dispersion_cutoff',       Set(60, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    # In the shifted case, it annoyingly takes one more iteration (12 vs 11)
    assertions = {PotentialType.TRUNCATED: OrderedDict([("energy", SetAssert(-19.399515443, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


def boron_nitride_hex(named_result='bn_hex') -> str:
    """
    Hexagonal.

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
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(6)),
                                                        ("energy", SetAssert(-9.426842352, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


def molybdenum_disulfide(named_result='MoS2') -> str:
    """
    Rhombohedral.

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
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


def gold_cadmium(named_result='CaAu') -> str:
    """
    Simple orthorhombic.
    TODO(Alex) Try to converge

    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + orthorhombic.simple_orthorhombic_cifs['gold_cadmium'].file
    crystal = cif_parser_wrapper(fname)
    arbitrary_shift = 0.18

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

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


def calcium_titanate(named_result='CaTiO3') -> str:
    """
    Simple orthorhombic.
    TODO(Alex) Try to converge

    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + orthorhombic.simple_orthorhombic_cifs['calcium_titanate'].file
    crystal = cif_parser_wrapper(fname)
    arbitrary_shift = 0.18

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

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


def copper_oxyfluoride(named_result='CuO2F') -> str:
    """
    Body-centred orthorhombic.

    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + orthorhombic.body_centred_orthorhombic_cifs['copper_oxyfluoride'].file
    crystal = cif_parser_wrapper(fname)
    arbitrary_shift = 0.38

    options = OrderedDict([
        ('h0_cutoff',               Set(60, 'bohr')),
        ('overlap_cutoff',          Set(60, 'bohr')),
        ('repulsive_cutoff',        Set(60, 'bohr')),
        ('dispersion_cutoff',       Set(60, 'bohr')),
        ('ewald_real_cutoff',       Set(60, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('symmetry_reduction',      Set(False)),
        ('solver',                  Set('SCC')),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(43)),
                                                        ("energy", SetAssert(-33.936868441, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))



def aluminium_titanate(named_result='TiAl2O5') -> str:
    """
    Base-centred C orthorhombic.
    # NOTE: Can't converge

    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + orthorhombic.base_centred_orthorhombic_cifs['aluminium_titanate'].file
    crystal = cif_parser_wrapper(fname)
    # Will require wrap_atoms
    arbitrary_shift = 0.38

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

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


def titanium_disilicide(named_result='TiSi2') -> str:
    """
    Face-centred orthorhombic.
    # NOTE: Can't converge

    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + orthorhombic.face_centred_orthorhombic_cifs['titanium_disilicide'].file
    crystal = cif_parser_wrapper(fname)
    # Will require wrap_atoms
    arbitrary_shift = 0.38

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

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


# TODO Find: simple monoclinic


# TODO Find: base-centred C monoclinic


def malh_perovskite(named_result='malh_perovskite'):
    """
    Second simple cubic, with several atoms in the central cell.
    # TODO Relaxed structure looks triclinic from lattice parameters
    # Space group suggests monoclinic
    # Need to modify alpha angle to 89.9 degrees but it does converge for truncated potential

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
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))


# TODO Add Triclinic


