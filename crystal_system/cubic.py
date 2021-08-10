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
cubic_cifs = {'alpha-N2': FileUrl("cifs/cubic/Simple/alpha-N2/N2_mp-25_primitive.cif",
                                  "https://materialsproject.org/materials/mp-25/"),
              'CH3NH3PbI3': FileUrl("cifs/cubic/Simple/2014_cubic_halides_PBEsol/CH3NH3PbI3.cif",
                                    "https://github.com/WMD-group/hybrid-perovskites/tree/master/")
              }

# Any space group beginning with I
bcc_cifs = {'potassium': FileUrl("cifs/cubic/BCC/potassium/K_mp-58_primitive.cif",
                                 "https://materialsproject.org/materials/mp-58/"),
            'sio2': FileUrl("cifs/cubic/BCC/sio2/SiO2_mp-1188220_primitive.cif",
                            "https://materialsproject.org/materials/mp-1188220/"),

            'nickel_triphosphide': FileUrl("cifs/cubic/BCC/nickel_triphosphide/NiP3_mp-2301_primitive.cif",
                        "https://materialsproject.org/materials/mp-2301")}

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
                                       "https://materialsproject.org/materials/mp-23703/"),
            'silicon': FileUrl("cifs/cubic/FCC/Si/Si_mp-149_primitive.cif",
                               "https://materialsproject.org/materials/mp-149/")
            }

# FCC but the lattice vectors are cubic
conventional_fcc_cifs = {
            'copper': FileUrl("cifs/cubic/FCC/Cu/Cu_mp-30_conventional_standard.cif",
                              "https://materialsproject.org/materials/mp-30/"),
            'sodium_chloride': FileUrl("cifs/cubic/FCC/NaCl/NaCl_mp-22862_conventional_standard.cif",
                                        "https://materialsproject.org/materials/mp-22862/"),
            'magnesium_oxide': FileUrl("cifs/cubic/FCC/MgO/MgO_mp-1265_conventional_standard.cif",
                                       "https://materialsproject.org/materials/mp-1265/"),
            'silicon': FileUrl("cifs/cubic/FCC/Si/Si_mp-149_conventional_standard.cif",
                               "https://materialsproject.org/materials/mp-149/")
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
            'bravais': bravais,
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