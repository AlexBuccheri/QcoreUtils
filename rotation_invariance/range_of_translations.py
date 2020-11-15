"""
Perform a range of translations for a given crystal.
Used BCC SiO2 as an example
One has to go crazy with cut-offs to get something close to a consistent
total energy (margin of 5.e-6 Ha) for the final shift of 0.25

This is demonstrated in the Gamma point calculation
"""

from collections import OrderedDict

from src.xtb_potential import xtb_potential_str, PotentialType
from src.utils import Set, SetAssert
from crystal_system import cubic, tetragonal, hexagonal, trigonal, orthorhombic
from src.pymatgen_wrappers import cif_parser_wrapper
from translation_invariance.string_generator import xtb_translational_invariance_string as translation_string


# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED


def sio2_zeolite(named_result='SiO2') -> str:
    """
    BCC structure
    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    arbitrary_shifts = [0., 0.1, 0.2, 0.24, 0.25]

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

    comments = ['! Fractional shift of '] * len(arbitrary_shifts)

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shifts, named_result,
                              comments=[x + str(arbitrary_shifts[i]) for i, x in enumerate(comments)])

def sio2_zeolite_gamma_point(named_result='SiO2') -> str:
    """
    BCC structure
    :param named_result: named result string
    :return: Input for testing translational invariance
    """
    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    arbitrary_shifts = [0., 0.1, 0.2, 0.24, 0.25]

    options = OrderedDict([
        ('h0_cutoff',               Set(60, 'bohr')),
        ('overlap_cutoff',          Set(60, 'bohr')),
        ('repulsive_cutoff',        Set(60, 'bohr')),
        ('dispersion_cutoff',       Set(60, 'bohr')),
        ('ewald_real_cutoff',       Set(60, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('potential_type',          Set(xtb_potential_str(potential_type))),
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(14)),
                                              ("energy", SetAssert(-69.074350228, 5.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    comments = ['! Fractional shift of '] * len(arbitrary_shifts)

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shifts, named_result,
                              comments=[x + str(arbitrary_shifts[i]) for i, x in enumerate(comments)])

print(sio2_zeolite_gamma_point())
quit()


def tio2_rutile_gamma(named_result='rutile') -> str:
    """
    Simple tetragonal
    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    crystal = tetragonal.tio2_rutile()
    arbitrary_shift = 0.24  #0.134

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(40, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("n_iter", SetAssert(12)),
                                                        ("energy", SetAssert(-22.578519079, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    return translation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result)