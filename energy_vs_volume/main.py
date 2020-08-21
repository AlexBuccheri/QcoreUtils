"""  Compute energy vs volume curves """
import numpy as np
import collections
import matplotlib.pyplot as plt


from crystal_system import cubic
from src.pymatgen_wrappers import cif_parser_wrapper
from src.utils import Set
from src import qcore_input_strings as qcore_input
from src.run_qcore import run_qcore

# def atomic_positions(crystal: dict, lattice_constant_factors: np.ndarray) -> tuple:
#     """
#     Get atomic positions for each choice of lattice constant
#
#     Expects simple cubic systems, although this would work in any case,
#     I'd just need to either a) suply lattice vectors or b) build lattice vectors
#     for any bravais lattice
#
#     :param crystal: dictionary containing crystal information from cif file
#     :param scaling: Lattice parameter scaling factors
#     :return: sets_of_positions : crystal positions for each lattice constant with units
#     of crystal['lattice_parameters']['a'].unit
#     """
#
#     assert len(crystal['lattice_parameters'].keys()) == 1, \
#         "Must only have one lattice parameter to uniformly increase volume"
#
#     keys = [key for key in crystal.keys()]
#     assert 'fractional' in keys, "Positions should be in fractional coordinates"
#
#     bulk_positions = np.asarray(crystal['fractional'])
#     a_bulk = [a.value for a in crystal['lattice_parameters'].values()][0]
#
#     sets_of_positions = []
#     for lattice_constant_factor in lattice_constant_factors:
#         lattice_constant = lattice_constant_factor * a_bulk
#         # Works because lattice is simple cubic, else matmul(lattice, bulk_positions)
#         positions = (lattice_constant * bulk_positions).tolist()
#         sets_of_positions.append(positions)
#
#     return sets_of_positions


def cubic_lattice_constants(crystal: dict, lattice_constant_factors: np.ndarray):
    """
    Get lattice constants for energy vs volume curves of cubic systems
    Performs checks and balances.
    Expects cubic systems i.e. one lattice parameter.

    :param crystal: dictionary containing crystal information from cif file
    :param scaling: Lattice parameter scaling factors
    :return:
    """
    assert crystal['bravais'] in ['cubic', 'bcc', 'fcc'], "crystal must be some form of cubic"

    assert len(crystal['lattice_parameters'].keys()) == 1, \
        "Must only have one lattice parameter to uniformly increase volume"

    keys = [key for key in crystal.keys()]
    assert 'fractional' in keys, "Positions should be in fractional coordinates"

    a_bulk = [a.value for a in crystal['lattice_parameters'].values()][0]
    lattice_constants = []
    for lattice_constant_factor in lattice_constant_factors:
        lattice_constants.append(lattice_constant_factor * a_bulk)

    return lattice_constants


def check_equal(lst) -> bool:
    return lst[1:] == lst[:-1]


def get_position_unit(crystal) -> str:
    units = []
    for parameter in crystal['lattice_parameters'].values():
        units.append(parameter.unit)
    assert check_equal(units), "Units of lattice parameters not consistent"
    return units[0]


def magnesium_oxide():
    file_name = '../' + cubic.conventional_fcc_cifs['magnesium_oxide'].file
    crystal = cif_parser_wrapper(file_name, fractional=True, is_primitive_cell=False, bravais='cubic')

    unit = get_position_unit(crystal)
    lattice_constant_factors = np.linspace(0.8, 1.2, 21, endpoint=True)
    lattice_constants = cubic_lattice_constants(crystal, lattice_constant_factors)

    named_result = 'mgo_volume'
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
        ('solver',                  Set('SCC'))
    ])

    total_energies = np.zeros(shape=(lattice_constant_factors.size))
    for i,al in enumerate(lattice_constants):
        lattice_factor = lattice_constant_factors[i]
        crystal['lattice_parameters']['a'] = Set(al, unit)
        input_string = qcore_input.xtb_input_string(crystal, settings, named_result=named_result)
        output = run_qcore(input_string)
        if not output:
            print('No result:', lattice_factor)
        else:
            total_energies[i] = output[named_result]['energy']
            # 4.19 ang is exp. See "Ab initio determination of the bulk properties of Mgo"
            print(lattice_factor, al, output[named_result]['energy'])

    return lattice_constant_factors, total_energies


lattice_constant_factors, total_energies = magnesium_oxide()
plt.plot(lattice_constant_factors, total_energies, label="MgO conventional")
plt.legend()
#plt.xlim(0, 1)
plt.xlabel("Lattice constant (bulk constant)")
plt.ylabel("Total energy (Ha)")
plt.show()