from src import qcore_input_strings as qcore_input, utils


def xtb_translational_invariance_string(
        crystal: dict, settings: dict, assertions: dict, shift: float, named_result: str, comments='') -> str:

    """
    Generate two input strings to test the translational invariance of periodix xTB in qCore.

    Parameters
    ----------
    crystal : dict
         crystal system: atomic positions (fractional or xyz), species, lattice parameters, space_group and n_atoms
    settings : dict
         xTB settings (excluding structure)
    assertions : dict
         assertions
    shift : float
        rigid shift to atomic positions, in fractional units
    named_result : str
        named result
    comments : str, optional
        File comments

    Results
    ----------
    qcore periodic xtb input string : str

    """

    # Named results
    nr_noshift = named_result + '_noshift'
    nr_shifted = named_result + '_shifted'

    # Lattice constants data should be in bohr to avoid a change in numbers
    # when qcore's internal 'angstrom_to_bohr' changes
    crystal['lattice_parameters'] = utils.angstrom_to_bohr(crystal['lattice_parameters'])

    # Original positions
    structure_no_shift = qcore_input.get_xtb_periodic_structure_string(crystal)

    # All other options
    other_options = qcore_input.commands_to_string(settings)

    # Assertions: Same values but different named results
    assertions_noshift = qcore_input.assertions_string(nr_noshift, assertions)
    assertions_shifted = qcore_input.assertions_string(nr_shifted, assertions)

    # Shifted positions
    crystal = utils.update_positions(crystal, shift)
    all_positions_in_cell = utils.check_fractional_positions(named_result, crystal['fractional'])
    assert all_positions_in_cell
    structure_shifted = qcore_input.get_xtb_periodic_structure_string(crystal)

    input1 = comments + "\n" + nr_noshift + ' := xtb(\n ' + structure_no_shift + '\n' + \
             other_options + '\n)\n' + assertions_noshift + '\n\n'
    input2 = comments + "\n" +  nr_shifted + ' := xtb(\n ' + structure_shifted + '\n' + \
             other_options + '\n)\n' + assertions_shifted + '\n\n'

    return input1 + input2