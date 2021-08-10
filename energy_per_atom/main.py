# Probably look into entos tests that use python
# and make this one that runs with ctest

from inspect import getmembers, isfunction

import energy_per_atom.crystals as crystals_module


def total_input_generator(crystals: list) -> str:
    """
    Retrieve functions from crystals_module with names
    that contain entries in crystals, then add input strings
    :param crystals: list of strings, containing crystal system names
    :return: total input string
    """
    functions = [item for item in getmembers(crystals_module) if isfunction(item[1])]
    total_string = ''
    for crystal in crystals:
        for (name, function) in functions:
            if crystal in name:
                total_string += function() + '\n'

    return total_string


crystals = ['silicon']

print(total_input_generator(crystals))





