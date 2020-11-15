import typing
import numpy as np
import math
import collections
import sys
import warnings

from src import unit_conversions

default_named_result = 'xtb_calc'

class Set():
    """ Container for value and unit  """
    def __init__(self, value: float, unit=None) -> None:
        self.value = value
        if unit is not None:
            self.unit = unit
        else:
            self.unit=''

    def update(self, value: float, unit=None):
        self.value = value
        if unit is not None:
            self.unit = unit


class FileUrl():
    """ Container for local file location and URL reference  """
    def __init__(self, file, url=None) -> None:
        self.file = file
        self.url = url
    # TODO(Alex) write this function
    # def get(self):
    #     if self.url is not None:
    #         # TODO(Alex) Download from the web
    #     else:
    #         exit("To get a crystal file from the web, a URL needs"
    #              "to be defined")
    #     return


class SetAssert:
    """
    Container for use with asserts, setting value and
    corresponding margin of error
    """
    def __init__(self, value, margin=''):
        self.value = value
        self.margin = margin


def angstrom_to_bohr(lattice_parameters: typing.Dict, precision=6) -> typing.Dict:

    """
    Convert lattice constant entries in lattice parameter dictionary
    to bohr, from angstrom.

    Parameters
    ----------
    lattice_parameters : dict
        Lattice constants and angles.
        Each dictionary value has a float value and str unit.
        Valid keys = a, b, c, alpha, beta, gamma
    precision : int, optional
        Number of decimal places to round constants to

    Returns
    -------
    Lattice_parameters : dict
        Lattice constants and angles.
        (a,b,c) will be in bohr.

    """

    ang_to_bohr = unit_conversions.angstrom_to_bohr
    for key, parameter in lattice_parameters.items():
        if parameter.unit == 'angstrom':
            new_value = round(parameter.value * ang_to_bohr, precision)
            lattice_parameters[key] = Set(new_value, 'bohr')
    return lattice_parameters


def list_to_string(mylist: list, joiner=',', precision=None) -> str:

    """
    Convert a list to string, with each entry separated by \p joiner.

    Parameters
    ----------
    mylist : list
       List to convert
    joiner : str, optional
       Character used to join entries of the list
    precision : int, optional
        Number of decimal places to round list entries to

    Returns
    -------
    mylist_str : str
        Elements of mylist in a string, with \p joiner
        between each entry.
    """

    # Don't round integers
    if all(isinstance(entry, int) for entry in mylist):
        assert precision is None

    if precision is None:
        return joiner.join([str(i) for i in mylist])
    else:
        return joiner.join([str(round(i, precision)) for i in mylist])


def ordered_dictionary(d: dict, order_by='key'):

    """
    Convert dictionary to ordered dictionary

    Parameters
    ----------
    d : dict
      dictionary
    order_by : str, optional
      Order to use for sorting

    Returns
    ----------
    Ordered dictionary : dict

    Notes
        Ref: https://docs.python.org/2/library/collections.html#collections.OrderedDict
    """

    # Ordered dictionaries supported by the standard
    python_version_minor = int(sys.version.split()[0][2])
    if python_version_minor >= 6:
        return d

    if order_by == 'key':
        index = 0
    elif order_by == 'value':
        index = 1
    else:
        exit("order_by only accepts 'key' and 'value'. Got ", order_by)

    return collections.OrderedDict(sorted(d.items(), key=lambda t: t[index]))


def generic_str(value: typing.Any) -> str:

    """
      Specific control over conversion of Any to sting.
      Largely superfluous as I want to keep brackets in my list strings

      Parameters
      ----------
      value : Any
        Value to convery

      Returns
      ----------
      string : str

      """

    if isinstance(value, str):
        return value
    # bool can also be evaluated as int, hence bool comes first
    elif isinstance(value, bool):
        return str(value).lower()
    elif isinstance(value, int) or isinstance(value, float):
        return str(value)
    elif isinstance(value, list):
        return str(value)
        # Want to keep brackets in
        #return list_to_string(value)
    else:
        quit("Have not implemented type ", type(value), "in generic_str")


def check_fractional_positions(named_result: str, positions: typing.List) -> bool:

    """
    Check all fractional coordinates are between 0 and 1.

    Parameters
    ----------
       named_result : str
          Identify what system the positions correspond to
       positions : list
          List of positions. Must be in fractional coordinates.

    Returns
    -------
        all_positions_in_cell : bool
            True => All positions in cell are within the limits of 0 and 1.
    """

    all_positions_in_cell = True

    for i,position in enumerate(positions):
        xyz = np.asarray(position)
        below_zero = np.where((xyz < 0))[0]
        exceed_one = np.where((xyz > 1))[0]
        if len(below_zero) > 0:
            print("Component/s of atom position " + str(i) + " of " + named_result + " below 0: " + str(xyz), file=sys.stderr)
            all_positions_in_cell = False
        elif len(exceed_one) > 0:
            print("Component/s of atom position " + str(i) + " of " + named_result + " exceed 1: " + str(xyz), file=sys.stderr)
            all_positions_in_cell = False

    return all_positions_in_cell


def get_positions_key(crystal: dict):

    """
    Get positions key from crystal dictionary

    Parameters
    ----------
    crystal : dict
       crystal system

    Returns
    -------
    positions_key : str
       Valid keys are 'fractional' and 'xyz'
    """

    if 'fractional' in crystal.keys() and 'xyz' in crystal.keys():
        exit('cannot have crystal coordinates stored as fractional and xyz')

    elif 'fractional' not in crystal.keys() and 'xyz' not in crystal.keys():
        exit('crystal position_key is missing. Must be fractional or xyz')

    return 'fractional' if 'fractional' in crystal.keys() else "xyz"

def substitute_positions_and_species(input_string: str, crystal: dict, position_key: str) -> str:

    """
    Fills in positions and species using string.format

    Assumes string substitutions are of the form {Xi:s} and {posi:s}
    for species and positions, respectively.

    Parameters
    ----------
    input_string : str
       String in which to make substitutions
    crystal : dict
       Crystal data
    position_key : str
       Units of atomic positions

    Returns
    -------
        Returns a string with the positions and species substituted in
    """

    assert position_key in ['fractional', 'xyz']

    strings = {}

    for i, position in enumerate(crystal[position_key]):
        strings['pos' + str(i)] = list_to_string(position, ', ', precision=5)

    for i, species in enumerate(crystal['species']):
        strings['X' + str(i)] = species

    # ** unpacks the dictionary
    return input_string.format(**strings)


def update_positions(crystal: dict, shift: float):

    """
    Applies a rigid shift to all positions of a crystal.


    Parameters
    ----------
    crystal : dict
       Crystal data
    shift : float
       Rigid shift, in units consistent with positions

    Returns
    -------
        crystal dictionary with updated atomic positions
    """

    position_key = get_positions_key(crystal)
    if position_key == 'fractional' and abs(shift) > 1:
        warnings.warn("magnitude of position shift exceeds |1| whilst "
                      "using fractional coordinates")

    new_positions = []
    for position in crystal[position_key]:
        new_positions.append([xyz + shift for xyz in position])
    crystal[position_key] = new_positions

    return crystal

def rotate_lattice_vector(crystal: dict, phi_z: float, phi_x: float, phi_y: float):

    """
    Applies a sequential rotation to the lattice vectors of a crystal.


    Parameters
    ----------
    crystal : dict
       Crystal data
    phi_z : float
       Perform rotation around the z axis
    phi_x : float
       Perform rotation around the x axis
    phi_y : float
       Perform rotation around the y axis

    Returns
    -------
        crystal dictionary with updated atomic positions
    """

    # switch from row-wise storage of vectors to column-wise
    new_positions = np.transpose(crystal["lattice_vectors"])

    # Rotate around z axis
    new_positions = np.dot(
                        np.array([[math.cos(phi_z), -math.sin(phi_z), 0],
                                  [math.sin(phi_z),  math.cos(phi_z), 0],
                                  [              0,                0, 1]]),
                        new_positions)

    # Rotate around x axis
    new_positions = np.dot(
                        np.array([[              1,                0,                0],
                                  [              0,  math.cos(phi_x), -math.sin(phi_x)],
                                  [              0,  math.sin(phi_x),  math.cos(phi_x)]]),
                        new_positions)

    #Rotate around y axis
    new_positions = np.dot(
                        np.array([[math.cos(phi_y),                0, -math.sin(phi_y)],
                                  [              0,                1,                0],
                                  [math.sin(phi_y),                0, math.cos(phi_y)]]),
                        new_positions)

    crystal["lattice_vectors"] = new_positions

    return crystal
