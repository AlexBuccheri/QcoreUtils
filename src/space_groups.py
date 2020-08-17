import numpy as np

# Useful references:
# https://en.wikipedia.org/wiki/List_of_space_groups
# https://en.wikipedia.org/wiki/Space_group#Table_of_space_groups_in_3_dimensions
# https://en.wikipedia.org/wiki/Crystal_system


crystal_families = ('triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'hexagonal', 'cubic')

crystal_systems = ('triclinic', 'monoclinic', 'orthorhombic', 'tetragonal', 'trigonal', 'hexagonal', 'cubic')

# In qcore, common crystal system abbreviations are used
qcore_bravais_lattices = ('triclinic',
                          'monoclinic',
                          'base_centred_monoclinic',
                          'orthorhombic',
                          'body_centred_orthorhombic',
                          'base_centred_orthorhombic',
                          'face_centred_orthorhombic',
                          'tetragonal',
                          'body_centred_tetragonal',
                          'rhombohedral',
                          'hexagonal',
                          'cubic',
                          'bcc',
                          'fcc')


# Note: I ASSUME A, B and C correspond to different base-centred configurations
# P = simple, is omitted in qcore inputs
configuration = {'P':'', 'I':'body_centred', 'F':'face_centred', 'A':'base_centred',
                 'B':'base_centred', 'C':'base_centred', 'R':'rhombohedral'}


def space_group_to_crystal_system(space_group_number: int) -> str:

    """
    For a given space group number, return the crystal system.

    Parameters
    ----------
    space_group_number : int
        space group number, ranging from 1 to 230

    Returns
    -------
    crystal_system : str
        String name for the crystal system

    """

    assert space_group_number > 0
    assert space_group_number <= 230

    triclinic_range    = np.arange(1,     2 + 1)
    monoclinic_range   = np.arange(3,    15 + 1)
    orthorhombic_range = np.arange(16,   74 + 1)
    tetragonal_range   = np.arange(75,  142 + 1)
    trigonal_range     = np.arange(143, 167 + 1)
    hexagonal_range    = np.arange(168, 194 + 1)
    cubic_range        = np.arange(195, 230 + 1)

    if space_group_number in triclinic_range:
        return 'triclinic'
    elif space_group_number in monoclinic_range:
        return 'monoclinic'
    elif space_group_number in orthorhombic_range:
        return 'orthorhombic'
    elif space_group_number in tetragonal_range:
        return 'tetragonal'
    elif space_group_number in trigonal_range:
        return 'trigonal'
    elif space_group_number in hexagonal_range:
        return 'hexagonal'
    elif space_group_number in cubic_range:
        return 'cubic'


# References:
# https://spglib.github.io/spglib/python-spglib.html#get-symmetry-from-database
# http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en
#
# #Generated with spglib:
#
# import spglib
#
# lattice_centring = {}
# for hall in range(1,530+1):
#     spacegroup_type = spglib.get_spacegroup_type(hall)
#     key = spacegroup_type['number']
#     value = spacegroup_type['international_short'][0]
#     lattice_centring[key] = value
#
# for key,value in lattice_centring.items():
#     print(str(key) + ":'" + str(value) + "',")
#
lattice_centring = { 1:'P', 2:'P', 3:'P', 4:'P', 5:'C', 6:'P', 7:'P', 8:'C', 9:'C', 10:'P',
                    11:'P', 12:'C', 13:'P', 14:'P', 15:'C', 16:'P', 17:'P', 18:'P', 19:'P', 20:'B',
                    21:'B', 22:'F', 23:'I', 24:'I', 25:'P', 26:'P', 27:'P', 28:'P', 29:'P', 30:'P',
                    31:'P', 32:'P', 33:'P', 34:'P', 35:'B', 36:'B', 37:'B', 38:'A', 39:'A', 40:'A',
                    41:'A', 42:'F', 43:'F', 44:'I', 45:'I', 46:'I', 47:'P', 48:'P', 49:'P', 50:'P',
                    51:'P', 52:'P', 53:'P', 54:'P', 55:'P', 56:'P', 57:'P', 58:'P', 59:'P', 60:'P',
                    61:'P', 62:'P', 63:'B', 64:'B', 65:'B', 66:'B', 67:'B', 68:'B', 69:'F', 70:'F',
                    71:'I', 72:'I', 73:'I', 74:'I', 75:'P', 76:'P', 77:'P', 78:'P', 79:'I', 80:'I',
                    81:'P', 82:'I', 83:'P', 84:'P', 85:'P', 86:'P', 87:'I', 88:'I', 89:'P', 90:'P',
                    91:'P', 92:'P', 93:'P', 94:'P', 95:'P', 96:'P', 97:'I', 98:'I', 99:'P', 100:'P',
                    101:'P', 102:'P', 103:'P', 104:'P', 105:'P', 106:'P', 107:'I', 108:'I', 109:'I', 110:'I',
                    111:'P', 112:'P', 113:'P', 114:'P', 115:'P', 116:'P', 117:'P', 118:'P', 119:'I', 120:'I',
                    121:'I', 122:'I', 123:'P', 124:'P', 125:'P', 126:'P', 127:'P', 128:'P', 129:'P', 130:'P',
                    131:'P', 132:'P', 133:'P', 134:'P', 135:'P', 136:'P', 137:'P', 138:'P', 139:'I', 140:'I',
                    141:'I', 142:'I', 143:'P', 144:'P', 145:'P', 146:'R', 147:'P', 148:'R', 149:'P', 150:'P',
                    151:'P', 152:'P', 153:'P', 154:'P', 155:'R', 156:'P', 157:'P', 158:'P', 159:'P', 160:'R',
                    161:'R', 162:'P', 163:'P', 164:'P', 165:'P', 166:'R', 167:'R', 168:'P', 169:'P', 170:'P',
                    171:'P', 172:'P', 173:'P', 174:'P', 175:'P', 176:'P', 177:'P', 178:'P', 179:'P', 180:'P',
                    181:'P', 182:'P', 183:'P', 184:'P', 185:'P', 186:'P', 187:'P', 188:'P', 189:'P', 190:'P',
                    191:'P', 192:'P', 193:'P', 194:'P', 195:'P', 196:'F', 197:'I', 198:'P', 199:'I', 200:'P',
                    201:'P', 202:'F', 203:'F', 204:'I', 205:'P', 206:'I', 207:'P', 208:'P', 209:'F', 210:'F',
                    211:'I', 212:'P', 213:'P', 214:'I', 215:'P', 216:'F', 217:'I', 218:'P', 219:'F', 220:'I',
                    221:'P', 222:'P', 223:'P', 224:'P', 225:'F', 226:'F', 227:'F', 228:'F', 229:'I', 230:'I'}


def space_group_to_bravais(space_group_number: int):

    """
    Get bravais lattice name from space group number

    :param space_group_number: int, space group number, ranging from 1 to 230
    :return: bravais: str, bravais lattice name, consistent with qcore input
    """

    assert space_group_number > 0
    assert space_group_number <= 230

    crystal_system = space_group_to_crystal_system(space_group_number)

    # Naming convention exceptions
    # qcore drops 'simple' prefix
    if lattice_centring[space_group_number] == 'P':
        bravais = crystal_system

    # Abbreviated names for cubic systems
    elif crystal_system == 'cubic':
        if lattice_centring[space_group_number] == 'P':
            bravais = 'cubic'
        if lattice_centring[space_group_number] == 'I':
            bravais = 'bcc'
        if lattice_centring[space_group_number] == 'F':
            bravais = 'fcc'

    # trigonal only in rhombohedral setting
    elif crystal_system == 'trigonal':
        if lattice_centring[space_group_number] == 'R':
            bravais = configuration[lattice_centring[space_group_number]]

    else:
        bravais = configuration[lattice_centring[space_group_number]] + "_" + crystal_system

    return bravais
