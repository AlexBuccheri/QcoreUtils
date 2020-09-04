"""
Get the CIF file of every elemental crystal and perform a QCore calc.

API Homepage:
https://materialsproject.org/infrastructure

API Examples:
https://gist.github.com/search?q=pymatgen+materials+api
"""

from pymatgen import MPRester
import pymatgen.io.cif

import os
import io

from src.pymatgen_wrappers import structure_parser_wrapper

# Atomic number to atomic symbol
an_to_symbol = {1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
                10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K',
                20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
                30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y',
                40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In',
                50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr',
                60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm',
                70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au',
                80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac',
                90: 'Th', 91: 'Pa', 92: 'U', 93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 99: 'Es',
                100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh', 108: 'Hs',
                109: 'Mt',
                110: 'Ds', 111: 'Rg', 112: 'Cn', 113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}


# Materials API key should be user-specific and not shared
# hence get from the env
def set_api_key():
    API_KEY = os.environ["MAPI_KEY"]
    if API_KEY is None:
        print("Each user needs to define their own API key as an env variable")
        quit("env variable 'MAPI_KEY' not found")
    return API_KEY


# Materials project supported properties
supported_properties = ("energy", "energy_per_atom", "volume",
                        "formation_energy_per_atom", "nsites",
                        "unit_cell_formula", "pretty_formula",
                        "is_hubbard", "elements", "nelements",
                        "e_above_hull", "hubbards", "is_compatible",
                        "spacegroup", "task_ids", "band_gap", "density",
                        "icsd_id", "icsd_ids", "cif", "total_magnetization",
                        "material_id", "oxide_type", "tags", "elasticity")

supported_task_properties = ("energy", "energy_per_atom", "volume",
                             "formation_energy_per_atom", "nsites",
                             "unit_cell_formula", "pretty_formula",
                             "is_hubbard",
                             "elements", "nelements", "e_above_hull",
                             "hubbards",
                             "is_compatible", "spacegroup",
                             "band_gap", "density", "icsd_id", "cif")


API_KEY = set_api_key()

n_elements = 1


# This does exactly what I want
# Will need to add something that goes through the entries and checks the number of atoms/space group in each
# as clearly a mix of elements
with MPRester(API_KEY) as mat_project:
    symbol = 'K'
    mp_computed = mat_project.get_entries(symbol, conventional_unit_cell=False, inc_structure="initial")[0]
    print(type(mp_computed.structure))
    print(mp_computed.structure)
    crystal = structure_parser_wrapper(mp_computed.structure, is_primitive_cell=True, fractional=True, remove_unused_parameters=True)
    print(crystal)


# def get_structure_from_mp(API_KEY: str, formula: str,
#                           is_primitive_cell=True,
#                           fractional=True,
#                           remove_unused_parameters=True):
#     """
#
#     materials_string = 'Cd-Se' for CdSe
#     Database syntax is MongoDB
#
#     See these great refs for scraping cifs from the database:
#     https://matsci.org/t/downloading-all-piezoelectric-materials/71/2
#     https://github.com/dwinston/mpbinder/blob/754f29537d417c7b67d6ee53ed630a2787b7149d/get_zintl_cifs.ipynb
#     """
#
#     # Convert formula to required material string i.e. CdSe to Cd-Se
#     material_string = formula[0]
#     for letter in formula[1:]:
#         if letter.isupper():
#             material_string += '-'
#         material_string += letter
#
#     mpr = MPRester(API_KEY)
#     docs = mpr.query({'chemsys': {'$in': [material_string]}}, ['material_id', 'pretty_formula', 'cif'])
#     print(len(docs))
#     cif_string = docs[0]['cif']
#
#     wrapper = io.TextIOWrapper(
#         io.BytesIO(),
#         encoding='cp1252',
#         line_buffering=True,
#     )
#     wrapper.write(cif_string)
#     parser = pymatgen.io.cif.CifParser(wrapper)
#
#     # parser = pymatgen.io.cif.CifParser()
#     # parser.from_string(cif_string)
#
#     print(parser.get_structures())
#     #structure = parser.get_structures(primitive=True)[0]
#     #print(structure)

# get_structure_from_mp(API_KEY, 'CdSe')

# for atomic_number in range(1, n_elements + 1):
#     with MPRester(API_KEY) as mat_project:
#         #an_to_symbol[atomic_number]
#         symbol = 'CdSe'
#         data = mat_project.get_task_data(symbol)[0]
#         print(data)
#         #print(mat_project.get_structures(symbol, final=False)[0])
#         mp_structure = mat_project.get_entries(symbol, conventional_unit_cell=False)[0]
#         print(mp_structure)
