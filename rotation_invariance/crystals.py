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
from rotation_invariance.string_generator \
    import xtb_rotational_invariance_string as rotation_string, \
    comments_list

import numpy as np


# Variable used by WHOLE script
potential_type = PotentialType.TRUNCATED

def sio2_zeolite(named_result='SiO2') -> str:
    """
    BCC structure
    Notes:
    If one exceeds 0.24, one has to greatly increase the cut-offs to
    get consistent total energies.
    See range_of_rotations.py for more exploration with this system

    :param named_result: named result string
    :return: Input for testing translational invariance
    """

    fname = '../' + cubic.bcc_cifs['sio2'].file
    crystal = cif_parser_wrapper(fname)
    assert crystal['bravais'] == 'bcc', "crystal not simple cubic"

    arbitrary_rotation = np.array([[1.3, 1.6, 1.2]])

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

    return rotation_string(crystal, options, assertions[potential_type],
                              arbitrary_rotation, named_result, comments=comments_list(arbitrary_rotation))

def tio2_rutile(named_result='rutile') -> str:
    """
    Simple tetragonal
    Notes:
    Similar issue to SiO2. If one goes to a larger translation where
    wrapping is required (0.24), then one has to turn up the cut-offs
    to obtain agreement to 7 d.p. in the total_energy.
    See the gamma_point input in range_of_/rotations.py

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

    return rotation_string(crystal, options, assertions[potential_type],
                              arbitrary_shift, named_result,
                              comments=comments_list(arbitrary_shift))



# TODO Find: simple monoclinic

# TODO Find: base-centred C monoclinic

# TODO Add Triclinic


