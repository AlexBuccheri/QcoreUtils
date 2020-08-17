"""
 Module containing trigonal crystal dictionaries:
  key : str
    crystal name
  value : str
    file path to cif
"""

from src.utils import Set

# Cubic crystals by bravais lattice
# Space groups: 143 - 167.

# Any space group beginning with P
hexagonal_cifs = {}

# Any space group beginning with R
rhomohedral_cifs = {}


# Crystals without cif files
# Avoid additional tabulating and use cif_parser_wrapper to generate
# these dictionaries.

def MoS2():

    """
    Molybdenum disulfide

    Notes
      Rhombohedral setting. Space group: 160
      Indirect band gap: 1.228 eV
      https://materialsproject.org/materials/mp-1434/
    """

    positions = [[0.000000, 0.000000, 0.000000],
                 [0.254176, 0.254176, 0.254176],
                 [0.412458, 0.412458, 0.412458]]
    species = ['Mo', 'S', 'S']
    bravais = 'rhombohedral'
    space_group = 160
    lattice_parameters = {'a': Set(6.856, 'angstrom'), 'alpha': Set(26.94, 'degree')}
    data = {'fractional': positions,
            'species': species,
            'lattice_parameters': lattice_parameters,
            'space_group': ('',space_group),
            'n_atoms': len(species)}
    return data
