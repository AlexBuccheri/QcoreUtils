"""
 Module containing orthrhombic crystal dictionaries with the signature:
  key : str
    crystal name
  value : str
    file path to cif

  calcium titanate, and aluminum hexathiohypodiphosphate details are tabulated.
"""

from src.utils import Set, FileUrl

# Orthorhombic crystals by bravais lattice
# Space groups: 16 - 74

# Any space group beginning with P
simple_orthorhombic_cifs = {'gold_cadmium': FileUrl("cifs/orthorhombic/Simple/CdAu/CdAu_mp-1404_primitive.cif",
                                                    "https://materialsproject.org/materials/mp-1404"),
                            'calcium_titanate': FileUrl("cifs/orthorhombic/Simple/CaTiO3/CaTiO3_mp-4019_primitive.cif"
                                                        "")
                            }

# Any space group beginning with C
# Need in C-configuration to be consistent with qCore lattice vectors
# and tabulated k-points
base_centred_orthorhombic_cifs = {'aluminium_titanate': FileUrl("cifs/orthorhombic/base_centred_C/TiAl2O5_mp-4930_primitive.cif",
                                                                "https://materialsproject.org/materials/mp-4930")}

# Any space group beginning with I
body_centred_orthorhombic_cifs = {'copper_oxyfluoride': FileUrl("cifs/orthorhombic/body_centred/CuO2F_mp-997107_primitive.cif",
                                                                "https://materialsproject.org/materials/mp-997107")
                                  }

face_centred_orthorhombic_cifs = {'titanium_disilicide': FileUrl("cifs/orthorhombic/face_centred/TiSi2_mp-2582_primitive.cif",
                                                                 "https://materialsproject.org/materials/mp-2582")}


# Crystals without cif files
# Avoid additional tabulating and use cif_parser_wrapper to generate
# these dictionaries.

def calcium_titanate():
    """
    calcium titanate

    Notes
      Ref: https://materialsproject.org/materials/mp-4019/
      Space group: Pnma [62]
    """

    positions = [[0.991521, 0.044799, 0.750000],
                 [0.491521, 0.455201, 0.250000],
                 [0.508479, 0.544799, 0.750000],
                 [0.008479, 0.955201, 0.250000],
                 [0.500000, 0.000000, 0.500000],
                 [0.000000, 0.500000, 0.500000],
                 [0.000000, 0.500000, 0.000000],
                 [0.500000, 0.000000, 0.000000],
                 [0.921935, 0.520580, 0.250000],
                 [0.421935, 0.979420, 0.750000],
                 [0.578065, 0.020580, 0.250000],
                 [0.078065, 0.479420, 0.750000],
                 [0.707456, 0.291917, 0.959281],
                 [0.207456, 0.208083, 0.040719],
                 [0.792544, 0.791917, 0.540719],
                 [0.292544, 0.708083, 0.459281],
                 [0.707456, 0.291917, 0.540719],
                 [0.207456, 0.208083, 0.459281],
                 [0.292544, 0.708083, 0.040719],
                 [0.792544, 0.791917, 0.959281]]

    species = ['Ca','Ca','Ca','Ca','Ti','Ti','Ti','Ti',
               'O ','O ','O ','O ','O ','O ','O ','O ','O ','O ','O ','O ']

    bravais = 'orthorhombic'

    space_group = 62
    lattice_parameters = {'a': Set(5.40444906, 'angstrom'),
                          'b': Set(5.51303112, 'angstrom'),
                          'c': Set(7.69713264, 'angstrom')}
    data = {'fractional': positions,
            'species': species,
            'lattice_parameters': lattice_parameters,
            'space_group': ('', space_group),
            'n_atoms': len(species)}

    return data


def aluminum_hexathiohypodiphosphate():
    """
    AlPS4.

    Notes
      Simple orthorhombic. Space group: P222 [16]
      Indirect band gap: 2.634 eV
      Ref: https://materialsproject.org/materials/mp-27462/
    """

    positions = [[0.000000, 0.000000, 0.000000],
                 [0.500000, 0.000000, 0.500000],
                 [0.000000, 0.500000, 0.000000],
                 [0.000000, 0.000000, 0.500000],
                 [0.197847, 0.276435, 0.101916],
                 [0.197847, 0.723565, 0.898084],
                 [0.802153, 0.276435, 0.898084],
                 [0.802153, 0.723565, 0.101916],
                 [0.776404, 0.800507, 0.601208],
                 [0.776404, 0.199493, 0.398792],
                 [0.223596, 0.800507, 0.398792],
                 [0.223596, 0.199493, 0.601208]]

    species = ['Al','Al','P','P','S','S','S','S','S','S','S','S']

    bravais = 'orthorhombic'

    space_group = 16
    lattice_parameters = {'a': Set(5.71230345, 'angstrom'),
                          'b': Set(5.71644625, 'angstrom'),
                          'c': Set(11.46678755,'angstrom')}
    data = {'fractional': positions,
            'species': species,
            'lattice_parameters': lattice_parameters,
            'space_group': ('', space_group),
            'n_atoms': len(species)}

    return data

