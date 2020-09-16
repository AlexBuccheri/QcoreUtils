"""
 Module containing cubic crystal dictionaries with the signature:
  key : str
    crystal name
  value : str
    file path to cif

  boron nitride is tabulated.
"""

from src.utils import Set

# Hexagonal crystals by bravais lattice
# Space groups: 168 - 194
# All space groups in this system begin with P

hexagonal_cifs = {}


# Crystals without cif files
# Avoid additional tabulating and use cif_parser_wrapper to generate
# these dictionaries.

def boron_nitride():

    """
    hexagonal boron nitride.

    Notes
      Space group: P63/mmc [194]
      Primitive lattice vectors and atomic basis
      Indirect and gap: 4.482 eV
      https://materialsproject.org/materials/mp-984/
    """

    positions = [[1/3, 2/3, 1/4],
                 [2/3, 1/4, 3/4],
                 [1/3, 2/3, 3/4],
                 [2/3, 1/3, 1/4]]
    species = ['B', 'B', 'N', 'N']
    bravais = 'hexagonal'
    space_group = 194
    lattice_parameters = {'a': Set(2.51242804, 'angstrom'), 'c': Set(7.70726501, 'angstrom')}
    data = {'fractional':positions,
            'species':species,
            'lattice_parameters':lattice_parameters,
            'space_group': ('', space_group),
            'bravais': bravais,
            'n_atoms': 4}
    return data


