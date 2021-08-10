""" Functions that generate partial input strings in qCore's
    parser language
"""

import typing

from src import space_groups, utils
from src.xtb_potential import PotentialType


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


def get_xtb_periodic_structure_string(crystal: dict,
                                      options=None) -> str:

    """
    For a periodic crystal, create atoms and lattice strings, fill
    in the atom species and positions, and concatenate the strings

    Parameters
    ----------
    crystal : dict
       Crystal data
    options : Optional, dict
       Options for structure command

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
    if options is None:
        return structure_string(atoms_str, lattice_str)
    else:
        options_str = option_to_string(options)
        return structure_string(atoms_str, lattice_str, options_str)


def option_to_string(options: dict, indent=1) -> str:

    """
    Create a string of qcore input options

    Parameters
    ----------
    options : dict
       A dictionary of qcore commands (as keys) and values with units (as value).
       For example, option, value = commands.items()
       and value is an object accessed as value.value and value.unit

    Returns
    -------
    commands_string : str
       String of input entries, each of the form: option = value unit \n

    """

    generic_str = utils.generic_str

    if isinstance(indent, int):
        indent = ' ' * indent
    else:
        assert isinstance(indent, (str, int)), "indent must be str or int"

    commands_string = ""
    for option, rhs in options.items():
        assert isinstance(option, str), "option isn't a string - passing option incorrectly"
        commands_string += indent + option + " = " + generic_str(rhs.value) + " " + rhs.unit + '\n'

    return commands_string


# Could probably use .join or something
def sum_strings(strings):
    summed = ''
    for string in strings:
        summed += string
    return summed


#TODO(Alex) Use this for structure and lattice, too
# Would make more sense to nest in the dict like one does with the commands
def commands_to_string(commands: dict) -> str:
    """
    Assumes nesting with first key of ordered dictionary defining
    the outer-most command and last key defining the inner-most command
    """
    command_str = ""
    indent = [' ']
    for command_key, options in commands.items():
        assert isinstance(command_key, str), "commands key isn't a string - passing commands incorrectly"
        assert isinstance(options, dict), "commands value isn't a dictionary - passing commands incorrectly"
        command_str += sum_strings(indent) + command_key + "(\n"
        indent.append(' ' * len(command_key + "("))
        command_str += option_to_string(options, sum_strings(indent))

    # Close command braces
    if len(commands) == 1:
        # Use an exception if no nested sub-commands
        command_str += sum_strings(indent[0:2]) + ") \n"
    else:
        for i in range(0, len(commands)):
            command_str += sum_strings(indent[0:i-1]) + ") \n"

    return command_str


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
    for variable_name, key in assertions.items():
        value = str(key.value)
        margin = str(key.margin)
        if margin != '':
            assert_string += " assert(load = " + named_result + \
                             " variable = " + variable_name + " value = " + value + " margin = " + margin + ")\n"
        else:
            assert_string += " assert(load = " + named_result + \
                             " variable = " + variable_name + " value = " + value + ")\n"
    return assert_string


def xtb_input_string(
        crystal: dict,
        settings:     typing.Optional[dict] = None,
        sub_commands: typing.Optional[dict] = None,
        named_result: typing.Optional[str] = None,
        assertions:   typing.Optional[str] = None,
        comments:     typing.Optional[str] = '') -> str:

    structure = get_xtb_periodic_structure_string(crystal)

    options = option_to_string(settings) if settings is not None else ''

    sub_commands = commands_to_string(sub_commands) if sub_commands is not None else ''

    named_result = named_result if named_result is not None else utils.default_named_result

    assertions = assertions_string(named_result, assertions) if assertions is not None else ''

    input_string = comments + "\n" + named_result + ' := xtb(\n ' + structure + '\n' + \
             options + sub_commands + '\n)\n' + assertions + '\n\n'

    return input_string


def xtb_input_string_cleaner(
        crystal:      typing.Optional[dict]=None,
        sub_commands: typing.Optional[dict]=None,
        options:      typing.Optional[dict]=None,
        assertions:   typing.Optional[dict]=None,
        named_result: typing.Optional[str]=None,
        comments='') -> str:
    assert crystal is not None or sub_commands is not None, \
        "Expect either crystal or sub_commands containing structure "

    if sub_commands:
        quit("Need to implement proper treatment of sub-commands")

    if assertions:
        assert not any(key in PotentialType for key in assertions.keys()), \
            "assertions dictionary should be passed in with a choice made for xtb PotentialType" \
            "i.e. assertions[PotentialType]"

    crystal['lattice_parameters'] = utils.angstrom_to_bohr(crystal['lattice_parameters'])

    structure_str = get_xtb_periodic_structure_string(crystal)

    options_str = option_to_string(options) if options is not None else ''

    sub_commands = commands_to_string(sub_commands) if sub_commands is not None else ''

    assertions_str = assertions_string(named_result, assertions) \
        if assertions is not None else ''

    named_result = named_result if named_result is not None else utils.default_named_result

    input = comments + "\n" + named_result + ' := xtb(\n ' + structure_str + '\n' + \
            options_str + '\n)\n' + assertions_str + '\n\n'

    return input