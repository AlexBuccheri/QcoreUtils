"""
 Module containing tetragonal crystal dictionaries with the signature:
  key : str
    crystal name
  value : str
    file path to cif

  rutile and anatase details are tabulated.
"""

from src.utils import Set

# Cubic crystals by bravais lattice
# Space groups: 75 - 142.

# Any space group beginning with P
simple_tetragonal_cifs = {}

# Any space group beginning with I
body_centred_tetragonal_cifs = {}


# Crystals without cif files
# Avoid additional tabulating and use cif_parser_wrapper to generate
# these dictionaries.

def tio2_rutile() -> dict:

    """
    TiO2 Rutile.

    Notes
      Space group P42/mnm (136)
      Direct band gap: 1.781 eV
      Ref: https://materialsproject.org/materials/mp-2657/
    """
    positions = [[0.000000, 0.000000, 0.000000],
                 [0.500000, 0.500000, 0.500000],
                 [0.695526, 0.695526, 0.000000],
                 [0.304474, 0.304474, 0.000000],
                 [0.195526, 0.804474, 0.500000],
                 [0.804474, 0.195526, 0.500000]]
    species = ['Ti', 'Ti', 'O', 'O', 'O', 'O']
    bravais = 'tetragonal'
    space_group = 136
    lattice_parameters = {'a': Set(4.6068, 'angstrom'), 'c': Set(2.9916, 'angstrom')}
    data = {'fractional': positions,
            'species': species,
            'lattice_parameters': lattice_parameters,
            'space_group': ('', space_group),
            'n_atoms': len(species)}

    return data


def tio2_anatase() -> dict:

    """
    TiO2 Anatase. Space group:  I4_1/amd (141)

    Notes
      Indirect band gap: 2.062 eV
      Ref: https://materialsproject.org/materials/mp-390/#

      SPGLib confirms that this is the primitive cell and is space group 141,
      however, it looks body-centred cubic rather than body-centred tetragonal
      from visual inspection of the lattice vectors and parameters.
    """
    positions = [[0.500000, 0.500000, 0.000000],
                 [0.250000, 0.750000, 0.500000],
                 [0.456413, 0.956413, 0.500000],
                 [0.706413, 0.706413, 0.000000],
                 [0.043587, 0.543587, 0.500000],
                 [0.293587, 0.293587, 0.000000]]
    species = ['Ti', 'Ti', 'O', 'O', 'O', 'O']
    bravais = 'body_centred_tetragonal'
    space_group = 141
    lattice_parameters = {'a': Set(5.55734663, 'angstrom'), 'c': Set(5.55734663, 'angstrom')}
    data = {'fractional': positions,
            'species': species,
            'lattice_parameters': lattice_parameters,
            'space_group': ('', space_group),
            'n_atoms': len(species)}

    return data

