"""
 Module containing cubic crystal dictionaries with the signature:
  key : str
    crystal name
  value : str
    file path to cif

  silicon, diamond and germanium details are tabulated.
"""

from src.utils import Set, FileUrl


# Cubic crystals by bravais lattice
# Space groups: 195-230.

# Any space group beginning with P
cubic_cifs = {}

# Any space group beginning with I
bcc_cifs = {}

# Any space group beginning with F
fcc_cifs = {
            'copper': FileUrl("cifs/cubic/FCC/Cu/Cu_mp-30_primitive.cif",
                               "https://materialsproject.org/materials/mp-30/"),
            'pbs': FileUrl("cifs/cubic/FCC/PbS/PbS_mp-21276_primitive.cif",
                            "https://materialsproject.org/materials/mp-21276/"),
            'palladium': FileUrl("cifs/cubic/FCC/Pa/Pa_mp-10740_primitive.cif",
                                  "https://materialsproject.org/materials/mp-10740/"),
            'boron_nitride': FileUrl("cifs/cubic/FCC/BN-cubic/BN_mp-1639_primitive.cif",
                                      "https://materialsproject.org/materials/mp-1639/"),
            'sodium_chloride': FileUrl("cifs/cubic/FCC/NaCl/NaCl_mp-22862_primitive.cif",
                                        "https://materialsproject.org/materials/mp-22862/"),
            'magnesium_oxide': FileUrl("cifs/cubic/FCC/MgO/MgO_mp-1265_primitive.cif",
                                       "https://materialsproject.org/materials/mp-1265/"),
            'lithium_hydride': FileUrl("cifs/cubic/FCC/LiH/LiH_mp-23703_primitive.cif",
                                       "https://materialsproject.org/materials/mp-23703/")
            }

# FCC but the lattice vectors are cubic
conventional_fcc_cifs = {
            'copper': FileUrl("cifs/cubic/FCC/Cu/Cu_mp-30_conventional_standard.cif",
                              "https://materialsproject.org/materials/mp-30/"),
            'sodium_chloride': FileUrl("cifs/cubic/FCC/NaCl/NaCl_mp-22862_conventional_standard.cif",
                                        "https://materialsproject.org/materials/mp-22862/"),
            'magnesium_oxide': FileUrl("cifs/cubic/FCC/MgO/MgO_mp-1265_conventional_standard.cif",
                                       "https://materialsproject.org/materials/mp-1265/")
            }

# Crystals without cif files
# Avoid additional tabulating and use cif_parser_wrapper to generate
# these dictionaries.

def silicon():

    """
    Silicon Crystal

    Notes
      Ref: https://doi.org/10.1103/PhysRevB.24.6121
      Experimental value in Table I

    """

    positions = [[0.00, 0.00, 0.00],[0.25, 0.25, 0.25]]
    species = ['Si', 'Si']
    bravais = 'fcc'
    space_group = 227
    lattice_parameters = {'a': Set(5.429, 'angstrom')}

    data = {'fractional': positions,
            'species': species,
            'lattice_parameters': lattice_parameters,
            'space_group': ('', space_group),
            'n_atoms': 2}

    return data


def diamond():
    """
    Diamond Crystal

    Notes
      Ref: https://doi.org/10.1103/PhysRevB.24.6121
      Experimental value in Table I

    """
    positions = [[0.00, 0.00, 0.00],[0.25, 0.25, 0.25]]
    species = ['C', 'C']
    bravais = 'fcc'
    space_group = 227
    lattice_parameters = {'a':Set(3.567, 'angstrom')}

    data = {'fractional': positions,
            'species': species,
            'lattice_parameters': lattice_parameters,
            'space_group': ('', space_group),
            'n_atoms': 2}

    return data


def germanium():
    """
    Germanium Crystal

    Notes
      Ref: https://doi.org/10.1103/PhysRevB.24.6121
      Experimental value in Table I

    """
    positions = [[0.00, 0.00, 0.00],[0.25, 0.25, 0.25]]
    species = ['Ge', 'Ge']
    bravais = 'fcc'
    space_group = 227
    lattice_parameters = {'a':Set(5.652, 'angstrom')}

    data = {'fractional': positions,
            'species': species,
            'lattice_parameters': lattice_parameters,
            'space_group': ('', space_group),
            'n_atoms': 2}

    return data