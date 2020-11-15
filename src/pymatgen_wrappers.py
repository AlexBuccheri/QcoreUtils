"""  Module providing a wrapper for pymatgen's cif parse, and functions that
     extract data from  pymatgen's Structure object for use in qCore
"""

import typing
import os
import pymatgen.io.cif
from pymatgen.core.structure import Structure

from src.utils import Set
from src import space_groups

# If Alex's set of modules is available
python_paths = os.environ['PYTHONPATH'].split(os.pathsep)
sub_dirs = [path.split('/')[-1] for path in python_paths]
if 'Python' in sub_dirs:
    import modules.electronic_structure.structure.bravais as lattices


def lattice_parameters_from_cif(structure, length_unit='angstrom', angle_unit='degree') \
        -> typing.Dict:

    """
    Get lattice constants and angles from pymatgen structure object
    and assign units to them.

    Parameters
    ----------
    structure : str
    length_unit : str, optional,
        Unit of lattice constants. Valid options = ['angstrom', 'bohr']
        cif file lattice constants are angstrom by default (verify)
    angle_unit : str : optional
        Unit of lattice angles. Valid options = ['degree', 'radian']
        cif file lattice angles are degrees by default (verify)

    Returns
    -------
    lattice_parameters_with_units : dict
        Valid keys = a, b, c, alpha, beta, gamma
        Each dictionary value is an object with a float value and str unit.

    """

    assert isinstance(structure, Structure)
    assert length_unit in ['angstrom', 'bohr']
    assert angle_unit in ['degree', 'radian']

    lattice_parameters = structure.as_dict()['lattice']
    lattice_parameters.pop('matrix')
    lattice_parameters.pop('volume')

    lattice_parameters_with_units = {}
    for key, value in lattice_parameters.items():
        if key in ['a', 'b', 'c']:
            unit = length_unit
        else:
            unit = angle_unit
        lattice_parameters_with_units[key] = Set(value, unit)

    return lattice_parameters_with_units


def remove_superflous_parameters(lattice_parameters: dict, bravais: str) -> dict:

    """
    Removes unneeded lattice parameter entries.
    Required as qCore will exit if superfluous lattice information is provided in input.

    Parameters
    ----------
    lattice_parameters : dict
       lattice constants and angles retrieved from cif file
    bravais : str
        bravais lattice

    Returns
    -------
    required_lattice_parameters : dict
        Lattice constants and angles required by qCore
    """

    required_lattice_parameters = {}

    if bravais in ['cubic', 'bcc', 'fcc']:
        required_lattice_parameters['a'] = lattice_parameters.pop('a')

    elif bravais in ['tetragonal', 'body_centred_tetragonal', 'hexagonal']:
        required_lattice_parameters['a'] = lattice_parameters.pop('a')
        required_lattice_parameters['c'] = lattice_parameters.pop('c')

    elif bravais == 'rhombohedral':
        required_lattice_parameters['a'] = lattice_parameters.pop('a')
        required_lattice_parameters['alpha'] = lattice_parameters.pop('alpha')

    elif bravais in ['orthorhombic', 'body_centred_orthorhombic',
                     'base_centred_orthorhombic', 'face_centred_orthorhombic']:
        required_lattice_parameters['a'] = lattice_parameters.pop('a')
        required_lattice_parameters['b'] = lattice_parameters.pop('b')
        required_lattice_parameters['c'] = lattice_parameters.pop('c')

    elif bravais in ['monoclinic', 'base_centred_monoclinic']:
        required_lattice_parameters['a'] = lattice_parameters.pop('a')
        required_lattice_parameters['b'] = lattice_parameters.pop('b')
        required_lattice_parameters['c'] = lattice_parameters.pop('c')
        required_lattice_parameters['alpha'] = lattice_parameters.pop('alpha')

    elif bravais == 'triclinic':
        required_lattice_parameters = lattice_parameters

    else:
        exit('choice of bravais lattice is not valid: ', bravais)

    return required_lattice_parameters



def cif_parser_wrapper(fname:str, is_primitive_cell=True, fractional=True, bravais=None,
                       remove_unused_parameters=True, supercell_coefficients=None,
                       lattice_vectors=None) -> typing.Dict:
    """
    Wrapper for pymatgen's cif parser.

    Extracts:
      * Atomic positions : np.array converted to list
      * Species labels : tuple(str) converted to list
      * lattice constants and angles : dict
      * space group : tuple(str, int)

    Parameters
    ----------
    fname : str
       Cif file name
    is_primitive_cell : bool, optional
       Have pymatgen convert lattice to primitive cell
       Returns primitive unit cell by default
    fractional : bool, optional
       Extract atomic basis positions in fractional else angstrom
       Uses fractional as default
    bravais : str, optional
       bravais lattice type. Required if unit cell is not primitive
    remove_unused_parameters : bool, optional
       Remove lattice parameters that are not required in the Qcore input
    supercell_coefficients : list of 3 integers, optional
       Integers to expand cell to supercell
       See structure.make_supercell() in https://pymatgen.org/usage.html

    Returns
    -------
    crystal_data : dict
        Dictionary holding atomic basis positions, species labels,
        lattice parameters and the space group.

    Notes
      Online documentation for pymatgen Structure object:
      https://pymatgen.org/pymatgen.core.structure.html?highlight=structure#module-pymatgen.core.structure

    """

    assert fname[-3:] == 'cif'

    if not is_primitive_cell and bravais is None:
        # For example, for conventional FCC bravais = simple cubic,
        # whereas primitive FCC bravais = FCC
        exit("when cif is read in as a conventional cell, the user must specify bravais lattice "
             "as the bravais lattice determined from the space group will not necessarily correspond to the "
             "lattice vectors")

    parser = pymatgen.io.cif.CifParser(fname)
    structure = parser.get_structures(primitive=is_primitive_cell)[0]

    if supercell_coefficients:
        assert len(supercell_coefficients) == 3, "supercell_coefficients should have length 3"
        structure.make_supercell(supercell_coefficients)

    if fractional:
        position_key = 'fractional'
        positions = structure.frac_coords.tolist()
    else:
        position_key = 'xyz'
        positions = structure.cart_coords.tolist()

    sg = structure.get_space_group_info()
    if bravais is None:
        bravais = space_groups.space_group_to_bravais(sg[1])

    lattice_parameters = lattice_parameters_from_cif(structure)

    # Compare tabulated lattice angles against those computed from lattice vectors
    # Make sure they agree
    #if 'Python' in sub_dirs:
    # TODO(Alex) need to modify this function to take kargs
    #    lattice_fn = lattices.get_lattice_vectors(bravais, lattice_parameters)

    # Required for qCore input
    if remove_unused_parameters:
        lattice_parameters = remove_superflous_parameters(lattice_parameters, bravais)
    species = [element_enum.value for element_enum in structure.species]
    assert len(species) == len(positions)

    crystal_data = {position_key: positions,
                   'species': species,
                   'lattice_parameters': lattice_parameters,
                   'space_group': sg,
                   'bravais': bravais,
                   'n_atoms': len(species)}

    if lattice_vectors is not None:
        crystal_data['lattice_vectors'] = lattice_vectors

    return crystal_data


def structure_parser_wrapper(structure:Structure,
                             is_primitive_cell=True,
                             fractional=True,
                             remove_unused_parameters=True) -> typing.Dict:
    """
    Same as above but expects pymatgen.core.structure.Structure object

    """

    if is_primitive_cell:
        structure = structure.get_primitive_structure()

    if fractional:
        position_key = 'fractional'
        positions = structure.frac_coords.tolist()
    else:
        position_key = 'xyz'
        positions = structure.cart_coords.tolist()

    sg = structure.get_space_group_info()
    bravais = space_groups.space_group_to_bravais(sg[1])
    lattice_parameters = lattice_parameters_from_cif(structure)

    # Required for qCore input
    if remove_unused_parameters:
        lattice_parameters = remove_superflous_parameters(lattice_parameters, bravais)
    species = [element_enum.value for element_enum in structure.species]
    assert len(species) == len(positions)

    crystal_data = {position_key: positions,
                    'species': species,
                    'lattice_parameters': lattice_parameters,
                    'space_group': sg,
                    'bravais': bravais,
                    'n_atoms': len(species)}

    return crystal_data


