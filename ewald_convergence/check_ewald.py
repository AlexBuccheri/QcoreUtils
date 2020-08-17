
import numpy as np

from crystal_system.tetragonal import tio2_rutile
from ewald_convergence.crystals import *  # Import all crystals
from src.run_qcore import run_qcore

from modules.electronic_structure.structure import atoms
from modules.electronic_structure.structure import bravais, supercell

# ewald = {"real": 40, "reciprocal": 20, "alpha": 0.2}
# named_result = "conventional_mgo"
# input_string = convention_magnesium_oxide(named_result, ewald=ewald)
# print(input_string)
# quit()
# result = run_qcore(input_string)
#
# print(result[named_result].keys())

# https://en.wikipedia.org/wiki/Hartree
ha_to_ev = 27.211386245989

#TODO(Alex) Would be useful to generalise
def super_cell_from_crystal(rutile: dict, cell_integers: list) -> dict:
    assert rutile['lattice_parameters']['a'].unit == 'angstrom'

    # Extract data
    a = rutile['lattice_parameters']['a'].value
    c = rutile['lattice_parameters']['c'].value
    lattice = bravais.simple_tetragonal(a, c)
    positions_angstrom = [np.matmul(lattice, pos) for pos in rutile['fractional']]
    species = rutile['species']

    # Compute supercell
    unit_cell = [atoms.Atom(species[ia], positions_angstrom[ia]) for ia in range(0, len(species))]
    super_cell = supercell.build_supercell(unit_cell, supercell.translation_vectors(lattice, cell_integers))
    assert len(super_cell) == np.prod(cell_integers) * len(unit_cell)
    lattice_scell = supercell.get_supercell_vector(lattice, cell_integers)
    #write.xyz('rutile_222.xyz', super_cell)

    # Repackage the information
    inv_lattice_scell = np.linalg.inv(lattice_scell)
    positions_frac = [np.matmul(inv_lattice_scell, atom.position) for atom in super_cell]
    species_sc = [atom.species for atom in super_cell]
    a = np.linalg.norm(lattice_scell[:, 0])
    c = np.linalg.norm(lattice_scell[:, 2])
    lattice_parameters = {'a': Set(a, 'angstrom'), 'c': Set(c, 'angstrom')}
    data = {'fractional': positions_frac,
            'species': species_sc,
            'lattice_parameters': lattice_parameters,
            'space_group': ('', 136),  #rutile['space_group']
            'n_atoms': len(species_sc)}

    return data


def rutile_input(named_result, cell_integers: list) -> str:

    rutile = tio2_rutile()
    rutile_supercell = super_cell_from_crystal(rutile, cell_integers)
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(40, 'bohr')), # Based on forcing xTB part to 0 at 30 bohr
        ('ewald_reciprocal_cutoff', Set(10)),
        ('ewald_alpha',             Set(0.5)),
        ('monkhorst_pack',          Set([2, 2, 2])),  # For speed
        ('symmetry_reduction',      Set(False)),
        ('temperature',             Set(0, 'kelvin')),
        ('solver',                  'SCC')
    ])

    input_string = qcore_input.xtb_input_string(rutile_supercell, settings, named_result=named_result)
    return input_string



# Cut-off of Potential. Internal 30 bohr cut-off and 1 bohr
atoms_in_primitive = 6
sets_of_cell_integers =  [[2, 2, 2], [3, 3, 3], [4, 4, 4]]
print("30 bohr, 1 bohr smoothing")
for cell_integers in sets_of_cell_integers:
    named_result = 'rutile'
    input_string = rutile_input(named_result, cell_integers)
    print(input_string)
    quit()
    json_output = run_qcore(input_string)
    print(json_output.keys())
    total_energy = json_output[named_result]["energy"] #* ha_to_ev
    energy_per_atom = total_energy / (atoms_in_primitive * np.prod(cell_integers))
    print("Total energy (Ha), energy per atom (Ha): ", total_energy, energy_per_atom)
    quit()