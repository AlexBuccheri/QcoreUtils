"""
 Module containing cubic crystal dictionaries with the signature:
  key : str
    crystal name
  value : str
    file path to cif

"""

# Monoclinic crystals by bravais lattice
# Space groups: 3-15.

# Any space group beginning with P
simple_monoclinic_cifs = {}

# Any space group beginning with C:
# Need in C-configuration to be consistent with qCore lattice vectors
# and tabulated k-points
base_centred_c_monoclinic_cifs = {}


