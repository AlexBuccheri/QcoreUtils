""" Functions that generate partial input strings in qCore's
    parser language
"""

import typing

import space_groups  # space_group_to_bravais, qcore_bravais_lattices
import utils         # substitute_positions


def atoms_string(n_atoms: int, position_key: str) -> str:

    """
    Generate qcore 'atoms input string', with placeholders for species and positions.
    To subsequently be processed with python's .replace()

    Parameters
    ----------
     n_atoms : int
        Number of atoms
     position_key : str
        Command key. Positions in 'xyz' or 'fractional'

    Returns
    -------
    atoms_string : str
        Returns fractional or xyz string of atomic species labels plus positions.

    """

    assert position_key in ['fractional', 'xyz']

    # First line
    atoms_string = position_key + "= [['{X0:s}', {pos0:s}],\n"
    indent = ' ' * (atoms_string.find('[') + 1)

    if n_atoms == 1:
        atoms_string = atoms_string.rstrip('\n')[:-1]
        return atoms_string + "]\n"

    if n_atoms > 2:
        for ia in range(1, n_atoms - 1):
            atoms_string += indent + "['{X" + str(ia) + ":s}', {pos" + str(ia) + ":s}],\n"

    # End line
    atoms_string += indent + "['{X" + str(n_atoms-1) + ":s}', {pos" + str(n_atoms-1) + ":s}]]\n"

    return atoms_string


def lattice_string(lattice_parameters: typing.Dict, bravais=None, space_group=None, precision=5) -> str:

    """
    Generate qcore lattice input string

    Parameters
    ----------
     lattice_parameters : Dict
         Dictionary of lattice constants and angles, where each value is an object
         with member data .value and .unit
     bravais : str, optional
        bravais lattice type (14 bravais lattice types)
     space_group : int, optional
        Space group number, 1-230. Bravais lattice information contained in space group
    precision : int, default = 5
        How many decimal places to round lattice constants/angles to

    Returns
    -------
    lattice_string : str
        Returns lattice string of lattice constants and angles, with units

    """

    if bravais is None:
        assert type(space_group) == tuple
        sg_number = space_group[1]
        assert type(sg_number) == int
        assert sg_number > 0
        assert sg_number < 231
        bravais = space_groups.space_group_to_bravais(sg_number)

    if space_group is None and bravais is None:
        exit("Need to pass bravais or space_group number to qcore 'lattice_string'")

    assert bravais in space_groups.qcore_bravais_lattices

    lattice_string = "lattice( \n"
    indent = " " * (lattice_string.find(" ") + 1)

    for key, rhs in lattice_parameters.items():
        value = str(round(rhs.value, precision))
        lattice_string += indent + key + " = " + value + " " + rhs.unit + "\n"

    lattice_string += indent + "bravais = " + bravais + "\n" + indent[:-1] + ")\n"

    return lattice_string


def structure_string(*args: str):

    """
    Put atom, lattice and any other subcommand strings into structure command string

    Parameters
    ----------
     *args : str
        Any number of subcommand strings

    Returns
    -------
    structure_string : str
        Returns Structure command string for qcore

    """

    structure_string = "structure( \n"
    indent = " " * (structure_string.find("(") + 1)
    for arg in args:
        lines = arg.split("\n")
        for line in lines:
            structure_string += indent + line + "\n"

    structure_string += indent[:-1] + ")\n"

    return structure_string


def get_xtb_periodic_structure_string(crystal: dict) -> str:

    """
    For a periodic crystal, create atoms and lattice strings, fill
    in the atom species and positions, and concatenate the strings

    Parameters
    ----------
    crystal : dict
       Crystal data

    Returns
    -------
    structure_string : str
        Returns Structure command string for qcore

    """
    position_key = utils.get_positions_key(crystal)
    empty_atoms_str = atoms_string(crystal['n_atoms'], position_key=position_key)
    atoms_str = utils.substitute_positions_and_species(empty_atoms_str, crystal, position_key=position_key)
    if 'bravais' in crystal:
        lattice_str = lattice_string(crystal['lattice_parameters'], bravais=crystal['bravais'])
    else:
        lattice_str = lattice_string(crystal['lattice_parameters'], space_group=crystal['space_group'])
    return structure_string(atoms_str, lattice_str)


def commands_to_string(commands: dict, indent=1) -> str:

    """
    Create a string of qcore input commands

    Parameters
    ----------
    commands : dict
       A dictionary of qcore commands (as keys) and values with units (as value).
       For example, command, value = commands.items()
       and value is an object accessed as value.value and value.unit

    Returns
    -------
    commands_string : str
       String of input entries, each of the form: command = value unit \n

    Notes
      Could expand to work with nested dictionaries
    """

    generic_str = utils.generic_str
    indent = ' ' * indent
    commands_string = ""

    for command, rhs in commands.items():
        assert isinstance(command, str)
        commands_string += indent + command + " = " + generic_str(rhs.value) + " " + rhs.unit + '\n'

    return commands_string


def assertions_string(named_result: str, assertions: dict) -> str:

    """
    Create a string of qcore assert commands

    Parameters
    ----------
    named_result : str
       qcore named result
    assertions : dict
       A dictionary of assertions: variable_name (key) and variable value (value)

    Returns
    -------
    assert_string : str
       String of assertions, each of the form:
       assert(load = named_result variable = variable_name value = variable_value)

    """

    assert_string = ''
    for variable_name, variable_value in assertions.items():
        assert_string += " assert(load = " + named_result + \
                         " variable = " + variable_name + " value = " + str(variable_value) + ")\n"
    return assert_string


def xtb_input_string(
        crystal: dict, settings: dict, named_result: str, assertions=None, comments='') -> str:

    structure = get_xtb_periodic_structure_string(crystal)

    other_options = commands_to_string(settings)

    if assertions is not None:
        assertions = assertions_string(named_result, assertions)
    else:
        assertions = ''

    input_string = comments + "\n" + named_result + ' := xtb(\n ' + structure + '\n' + \
             other_options + '\n)\n' + assertions + '\n\n'

    return input_string