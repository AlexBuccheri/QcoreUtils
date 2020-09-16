from collections import OrderedDict
import copy

from src import qcore_input_strings as qcore_input, utils


def comments_list(shifts) -> list:
    """
    Create comments list containing comment strings for each choice of shift
    :param shifts: float or list of floats, shift(s) in fractional coordinates
    :return: comments list
    """
    assert isinstance(shifts, (float, list)), "shift must a float or list of floats"

    if isinstance(shifts, float):
        return ['! No shift', '! Fractional shift of ' + str(shifts)]

    elif isinstance(shifts, list):
        comments = []
        for shift in shifts:
            if shift == 0.:
                comments.append('! No shift')
            else:
                comments.append('! Fractional shift of ' + str(shifts))
        return comments


def xtb_input_string(
        crystal:      dict,
        options:      dict,
        assertions:   dict,
        named_result: str,
        sub_commands=None,
        comments=None) -> str:
    """
    Generate two input strings to test the translational invariance of periodic xTB in qCore.

    Parameters
    ----------
    crystal : dict
         crystal system: atomic positions (fractional or xyz), species, lattice parameters, space_group and n_atoms
    sub_commands : dict
         sub-commands excluding structure
    options : dict
         xTB options (excluding sub-commands, like structure)
    assertions : dict
         assertions
    shift : float
        rigid shift to atomic positions, in fractional units
    named_result : str
        named result

    Results
    ----------
    qcore periodic xtb input string : str

    """
    crystal['lattice_parameters'] = utils.angstrom_to_bohr(crystal['lattice_parameters'])
    all_positions_in_cell = utils.check_fractional_positions(named_result, crystal['fractional'])

    if all_positions_in_cell:
        structure_str = qcore_input.get_xtb_periodic_structure_string(crystal)
    else:
        structure_options = OrderedDict([('wrap_atoms', utils.Set(True))])
        structure_str = qcore_input.get_xtb_periodic_structure_string(crystal, structure_options)

    sub_commands_str = qcore_input.commands_to_string(sub_commands) if sub_commands is not None else ''
    options_str = qcore_input.option_to_string(options)
    assertions_str = qcore_input.assertions_string(named_result, assertions)
    comments_str = comments if comments else ''

    return comments_str + '\n' + named_result + ' := xtb(\n ' + structure_str + '\n' + \
             options_str + sub_commands_str + '\n)\n' + assertions_str + '\n'


def xtb_translational_invariance_string(
        crystal: dict,
        options: dict,
        assertions: dict,
        shift,
        named_result: str,
        comments=None) -> str:
    """
    Translational invariance of the total energy means that the assertions can be the same
    """
    assert isinstance(shift, (float, list)), "shift must a float or list of floats"

    if isinstance(shift, float):
        if comments:
            assert len(comments) == 2, "len(comments) != 2"
        else:
            comments = [''] * 2

        input_no_shift = xtb_input_string(crystal, options, assertions, named_result + '_no_shift', comments=comments[0])
        crystal = utils.update_positions(crystal, shift)
        input_shift = xtb_input_string(crystal, options, assertions, named_result + '_shift', comments=comments[1])
        return input_no_shift + '\n' + input_shift

    elif isinstance(shift, list):
        if comments:
            assert len(comments) == len(shift), "len(comments) != len(shift)"
        else:
            comments = [''] * len(shift)

        input = ''
        for i,s in enumerate(shift):
            assert isinstance(s, float), "fractional shift must a float or list of floats"
            updated_crystal = utils.update_positions(copy.deepcopy(crystal), s)
            input += xtb_input_string(updated_crystal, options, assertions,
                                      named_result + '_shift' + str(i), comments=comments[i]) + '\n'
        return input

