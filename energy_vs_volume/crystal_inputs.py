import numpy as np
import collections

from crystal_system import cubic
from src.pymatgen_wrappers import cif_parser_wrapper
from src.utils import Set
from src import qcore_input_strings as qcore_input
from src.run_qcore import run_qcore

from energy_vs_volume import ev_functions


# MgO, NaCl, Si, Diamond (carbon), Ge (although could be an issue)
# Add in a few more

def magnesium_oxide():
    file_name = '../' + cubic.conventional_fcc_cifs['magnesium_oxide'].file
    crystal = cif_parser_wrapper(file_name, fractional=True, is_primitive_cell=False, bravais='cubic')

    unit = ev_functions.get_position_unit(crystal)
    lattice_constant_factors = np.linspace(0.8, 1.2, 21, endpoint=True)
    lattice_constants = ev_functions.cubic_lattice_constants(crystal, lattice_constant_factors)

    named_result = 'mgo_volume'
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        # Ewald setting for hard-cutoff of potential at 30 bohr
        ('ewald_real_cutoff',       Set(40, 'bohr')),
        # Converged w.r.t. real-space value
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
