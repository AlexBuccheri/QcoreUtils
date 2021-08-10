"""  Compute energy vs volume curves """

import numpy as np
import matplotlib.pyplot as plt
import itertools
import scipy.optimize

def cubic_lattice_constants(crystal: dict, lattice_constant_factors: np.ndarray)\
        -> np.ndarray:
    """
    Get lattice constants for energy vs volume curves of cubic systems
    Performs checks and balances.
    Expects cubic systems i.e. one lattice parameter.

    :param crystal : dictionary containing crystal information from cif file
    :param scaling : Lattice parameter scaling factors
    :return : Array of lattice constants in whatever unit the crystal['lattice_parameters']
    are in (i.e. the unit of the cif file, which I believe is always angstrom)
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


def plot_energy_vs_volume(crystal_input, labels):
    lattice_constant_factors, total_energies = crystal_input()
    plt.plot(lattice_constant_factors, total_energies, label=labels['crystal_cell'])
    plt.legend()
    plt.xlabel("Lattice constant (bulk constant)")
    plt.ylabel("Total energy (Ha)")
    plt.show()

# TODO(Alex)
# Have something that quantifies the error too
# Equilibrium lattice
# Error in curve (delta)



def birch_murnaghan_relation(V:np.ndarray, E0: float, V0: float, B0: float, B1: float):
    """
    Third-order Birch-Murnaghan relation, relating energy to volume.

    Given by equation 2 in 10.1080/10408436.2013.772503 and references therein.
    See https://doi.org/10.3390/min9120745 for a thorough derivation.

    :param V: Volume data
    :param E0: Parameter. Energy at Equilibrium volume
    :param E0: Parameter. Equilibrium volume
    :param B0: Parameter. bulk modulus
    :param B1: Parameter. Derivative of the bulk modulus w.r.t. pressure
    :return: third-order Birch-Murnaghan relation
    """
    term1 = ((V0 / V)**(2/3) - 1)
    return E0 + (9 / 16 * V0 * B0) * ( (term1**3 * B1) + (term1**2 * (6 - 4 * (V0 / V)**(2/3))) )


def check_keys(a: dict, expected_keys: list):
    """

    :return:
    """
    assert type(a) == dict
    keys = [key for key in a.keys()]
    for key in keys:
        assert key in expected_keys


def fit_equation_of_state(data:dict, volume_grid:np.ndarray):
    """
    Given a set of volumes, return corresponding energies via
    fitting to the birch murnaghan relation

    Alternatively could have interpolated but birch murnaghan relation
    should provide an excellent fit to the data:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html#scipy.interpolate.interp1d

    :return: Energies corresponding to volumes in volume_grid, fitting parameters
    and their standard deviations
    """
    check_keys(data, ['E', 'V'])
    param_guesses = [np.amin(data['E']), np.amin(data['V']), 1., 1.]
    params, pcov = scipy.optimize.curve_fit(
        birch_murnaghan_relation, data['V'], data['E'], param_guesses)
    p_standard_deviations = np.sqrt(np.diag(pcov))
    return birch_murnaghan_relation(volume_grid, *params), params, p_standard_deviations


def equation_of_states_rms(data: dict, data_ref: dict):
    """
    Root mean squared error for two equations of state

    :param data: energy vs volume
    :param data_ref: reference energy vs volume
    :return: RMS error
    """
    check_keys(data, ['E', 'V'])
    check_keys(data_ref, ['E', 'V'])

    # evaluate both continuous E(V) functions on a regular grid
    # within a range over which data points exist
    V_min = max(np.amin(data['V']), np.amin(data_ref['V']))
    V_max = min(np.amax(data['V']), np.amax(data_ref['V']))
    volume_grid = np.linspace(V_min, V_max,  num=100, endpoint=True)

    e_computed, p, sd = fit_equation_of_state(data, volume_grid)
    e_reference, p_ref, sd_ref = fit_equation_of_state(data_ref, volume_grid)
    error()

    # TODO(Alex) Check parameter SDs
    return np.sqrt(np.mean((e_computed - e_reference) ** 2))



# def delta_error():
#     """
#
#     $ \Delta = $
#
#     The root mean squared (rms) energy difference between two E vs V curves
#     averaged over all elemental crystals. (i.e. over all􏰀 volumes).
#
#     Provides an intuitive measure of the energy distance between equations of state.
#
#     :return:
#     """


def i_plus_j_equal_to_n_plus_2(n_values, i_plus_j):
    """

    Reference:
    https://stackoverflow.com/questions/34517540/find-all-combinations-of-a-list-of-numbers-with-a-given-sum

    :param n_values:
    :param i_plus_j:
    :return:
    """
    i_j_values_equal_to_n_plus_2 = {}
    for n_plus_2 in n_values:
        # Might be better to not have this all in-place (i.e creating extra work)
        result = [seq for index in range(len(i_plus_j), 0, -1)
                  for seq in itertools.combinations(i_plus_j, index)
                  if sum(seq) == n_plus_2]
        # where list(set( )) removes duplicates
        meaningful_results = list(set([ij for ij in result if len(ij) == 2]))
        i_j_values_equal_to_n_plus_2[n_plus_2 - 2] = meaningful_results
    return i_j_values_equal_to_n_plus_2


def delta_coefficients(a, a_ref):
    """
    Calculate the coefficients:
      x_n = - \frac{3}{2n+1} \sum_{i + j = n + 2} (a_i - a_i^{ref}) (a_j - a_j^{ref})

    for use in the power series expansion of F(V), used in the evaluation of the
    Delta function, quantifying RMS error in two equations of state
    (energy vs volume curves)

    Error Estimates for Solid-State Density-Functional Theory Predictions:
    An Overview by Means of the Ground-State Elemental Crystals, Critical Reviews in
    Solid State and Materials Sciences, 39:1, 1-24, DOI: 10.1080/10408436.2013.772503

    :param a: coefficients of the Birch-Murnaghan equation in its polynomial form
              for the computed data
    :param a_ref: coefficients of the Birch-Murnaghan equation, for the reference data
    :return: Delta coefficients, x_n
    """
    assert len(a) == 4, "Expect coefficients [a_0, a_1, a_2, a_3, a_4]"

    i = (0, 1, 2, 3)
    j = (0, 1, 2, 3)
    i_plus_j = i + j
    n_values = (-2, -1, 0, 1, 2, 3, 4)
    n_plus_two = [n + 2 for n in n_values]
    # Could just compute once and tabulate
    ij_values = i_plus_j_equal_to_n_plus_2(n_plus_two, i_plus_j)

    x = {}
    for n, ij_permutations in ij_values.items():
        prefactor = - 3 / (2 * n + 1)
        x_value = 0
        for ij in ij_permutations:
            i, j = ij
            x_value += (a[i] - a_ref[i]) * (a[j] - a_ref[j])
        x[n] = prefactor * x_value

    return x










