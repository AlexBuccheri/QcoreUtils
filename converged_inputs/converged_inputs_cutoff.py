"""
Converged inputs using truncated real-space sum in periodic xTB

---------------
Semiconductors
---------------
IV:     Si, Ge, Diamond, alpha-Sn. All FCC. No charge transfer between atoms
III-V:  Boron nitride (cubic and hex), GaAs, GaP, InAs, InSb
II-VI:  CdSe (hex), CdS, CdTe, ZnS, ZnSe, ZnTe,
IV-VI:  PbS (all have comparably small band gaps)

-------------------------
Wide bandgap/insulators
-------------------------
Oxide: TiO2 rutile and anatase, ZnO
Insulators: NaCl, MgO

Layered: MoS2, Graphite buckled and unbuckled.
https://pubs.acs.org/doi/10.1021/acs.jpclett.0c00721

---------------
Metals
---------------
Metals should largely get covered via running elemental crystals.

---------------
Applications
---------------
Thin film: Cu2ZnSnS4 (CZTS), Cu2SnS3 (CTS), CdZnTe (CTZ)
Zeolite: SiO2
Perovskites: CH3NH3PbI3  Be careful of correlation (DFT+U), SO coupling, etc. May not be viable

-------------------
Elemental Crystals
-------------------
Can run with excessively high cut-offs, then increase them by 50% and show convergence (all automated)
Will be hard to automatically process/compare the results.

"""

import collections

from crystal_system import cubic, tetragonal, trigonal, orthorhombic, monoclinic
from src.pymatgen_wrappers import cif_parser_wrapper
from src import qcore_input_strings as qcore_input
from src.utils import Set, default_named_result
from src.run_qcore import run_qcore


# Simple cubic
# Alpha N2, potentially LaAlO3

#TODO(Alex) Converge
def alpha_N2(named_result='alpha_N2'):
    """
    Simple cubic example
    No charge transfer => Trivial to converge
    """
    fname = '../' + cubic.cubic_cifs['alpha-N2'].file
    crystal = cif_parser_wrapper(fname)
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(40, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(10)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    xtb_potential = collections.OrderedDict([
        ('potential_type', Set('truncated')),
        ('smoothing_range', Set(1, 'bohr'))
    ])

    sub_commands = collections.OrderedDict([
        ('xtb_potential', xtb_potential)
    ])

    input_string = qcore_input.xtb_input_string(crystal,
                                                settings=settings,
                                                sub_commands=sub_commands,
                                                named_result=named_result)
    return input_string



# BCC
# potassium (a metal), SiO2 zeolite

#TODO(Alex) Need to converge Ewald and MP
def primitive_potassium(named_result='primitive_potassium'):
    """
    BCC Example

    Notes:
    Converges but with this warning
    warning: solve(): system seems singular (rcond: 4.98607e-17);
    attempting approx solution
    13         7.217892848     -5.61e-12      5.55e-17        0.13

    SCF takes ~ 35 iterations => SCC faster
    symmetry reduction works

    :param named_result: named result
    :return:
    """
    fname = '../' + cubic.bcc_cifs['potassium'].file
    crystal = cif_parser_wrapper(fname)

    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(40, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(10)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin')),
        ('solver',                     Set('SCC'))
    ])

    xtb_potential = collections.OrderedDict([
        ('potential_type', Set('truncated')),
        ('smoothing_range', Set(1, 'bohr'))
    ])

    sub_commands = collections.OrderedDict([
        ('xtb_potential', xtb_potential)
    ])

    input_string = qcore_input.xtb_input_string(crystal,
                                                settings=settings,
                                                sub_commands=sub_commands,
                                                named_result=named_result)
    return input_string









# FCC
# si, ge, diamond
# sodium_chloride, magnesium_oxide, PbS, CdS





# Simple tetragonal
# rutile


# Body-centred tetragonal
#TODO(Alex) Need to converge Ewald and MP
def primitive_anatase(named_result='primitive_anatase'):
    """
    Anatase primitive unit cell

    SCC takes 26 iterations vs 11 SCF

    :param named_result: named result
    :return:
    """
    crystal = tetragonal.tio2_anatase()

    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(40, 'bohr')),
        ('ewald_reciprocal_cutoff', Set(10)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    xtb_potential = collections.OrderedDict([
        ('potential_type', Set('truncated')),
        ('smoothing_range', Set(1, 'bohr'))
    ])

    sub_commands = collections.OrderedDict([
        ('xtb_potential', xtb_potential)
    ])

    input_string = qcore_input.xtb_input_string(crystal,
                                                settings=settings,
                                                sub_commands=sub_commands,
                                                named_result=named_result)
    return input_string

print(primitive_anatase())


# hexagonal
# CdSe, BN, graphite


# rhombohedral


# simple ortho


# base_centred_ortho


# body_centred_ortho


# face_centred_ortho


# simple monoclinic


# base-centred monoclinic


# triclinic