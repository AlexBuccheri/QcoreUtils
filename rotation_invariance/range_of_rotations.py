"""
Perform Rotation of lattice vectors for a given crystal.
Used BCC SiO2 as an example

This is demonstrated in the Gamma point calculation
"""

from collections import OrderedDict

from src.xtb_potential import xtb_potential_str, PotentialType
from src.utils import Set, SetAssert
from crystal_system import cubic, tetragonal, hexagonal, trigonal, orthorhombic
from src.pymatgen_wrappers import cif_parser_wrapper
from rotation_invariance.string_generator import xtb_rotational_invariance_string as rotation_string

import math
import numpy as np

# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED


def sio2_zeolite(named_result='SiO2') -> str:
    """
    BCC structure
    :param named_result: named result string
    :return: Input for testing rotational invariance
    """
    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname,
                                 lattice_vectors=np.array([[-7.35759, 7.35759, 7.35759],
                                                           [ 7.35759,-7.35759, 7.35759],
                                                           [ 7.35759, 7.35759,-7.35759]]))
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    arbitrary_rotations = np.array([[0.0, 0.0, 0.0],
                                    [math.pi / 2.0, 0.0, 0.0],
                                    [0.0001, 0.0001, 0.0001],
                                    [math.pi / 6.0, math.pi / 3.0, math.pi / 4.0],
                                    [2 * math.pi, 2 * math.pi, 2 * math.pi],
                                    [1.1, 2.2, 0.7],
                                    [0.6, 0.25, 1.4]
                                   ])

    options = OrderedDict([
        ('h0_cutoff',               Set(80, 'bohr')),
        ('overlap_cutoff',          Set(80, 'bohr')),
        ('repulsive_cutoff',        Set(80, 'bohr')),
        ('ewald_real_cutoff',       Set(80, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('potential_type',          Set(xtb_potential_str(potential_type))),
    ])

    assertions = {
        PotentialType.TRUNCATED: OrderedDict([("energy", SetAssert(-69.056517465, 1.e-5))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    comments = ['! Rotation of '] * len(arbitrary_rotations)

    return rotation_string(crystal, options, assertions[potential_type],
                              arbitrary_rotations, named_result,
                              comments=[x + str(arbitrary_rotations[i]) for i, x in enumerate(comments)])

def sio2_zeolite_gamma_point(named_result='SiO2') -> str:
    """
    BCC structure
    :param named_result: named result string
    :return: Input for testing rotational invariance
    """
    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname,
                                 lattice_vectors=np.array([[-7.35759, 7.35759, 7.35759],
                                                           [ 7.35759,-7.35759, 7.35759],
                                                           [ 7.35759, 7.35759,-7.35759]]))
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    arbitrary_rotations = np.array([[0.0, 0.0, 0.0],
                                    [math.pi / 2.0, 0.0, 0.0],
                                    [0.0001, 0.0001, 0.0001],
                                    [math.pi / 6.0, math.pi / 3.0, math.pi / 4.0],
                                    [2 * math.pi, 2 * math.pi, 2 * math.pi],
                                    [1.1, 2.2, 0.7],
                                    [0.6, 0.25, 1.4]
                                   ])

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
        PotentialType.TRUNCATED: OrderedDict([("energy", SetAssert(-69.074350228, 5.e-6))]),
        PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                              ("energy", SetAssert(0, 0))])
    }

    comments = ['! Rotation of '] * len(arbitrary_rotations)

    return rotation_string(crystal, options, assertions[potential_type],
                              arbitrary_rotations, named_result,
                              comments=[x + str(arbitrary_rotations[i]) for i, x in enumerate(comments)])


def tio2_rutile_gamma(named_result='rutile') -> str:
    """
    Simple tetragonal
    :param named_result: named result string
    :return: Input for testing rotational invariance
    """

    crystal = tetragonal.tio2_rutile(lattice_vectors=8.70559 * np.eye(3))

    arbitrary_rotations = np.array([[0.0, 0.0, 0.0],
                                    [math.pi / 2.0, 0.0, 0.0],
                                    [0.0001, 0.0001, 0.0001],
                                    [math.pi / 6.0, math.pi / 3.0, math.pi / 4.0],
                                    [2 * math.pi, 2 * math.pi, 2 * math.pi],
                                    [1.1, 2.2, 0.7],
                                    [0.6, 0.25, 1.4]
                                   ])

    options = OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(40, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(1)),
        ('ewald_alpha',             Set(0.5)),
        ('potential_type',          Set(xtb_potential_str(potential_type)))
    ])

    assertions = {PotentialType.TRUNCATED: OrderedDict([("energy", SetAssert(-22.532483264, 1.e-6))]),
                  PotentialType.FULL:      OrderedDict([("n_iter", SetAssert(0, 0)),
                                                        ("energy", SetAssert(0, 0))])
                  }

    comments = ['! Rotation of '] * len(arbitrary_rotations)

    return rotation_string(crystal, options, assertions[potential_type],
                              arbitrary_rotations, named_result,
                              comments=[x + str(arbitrary_rotations[i]) for i, x in enumerate(comments)])
