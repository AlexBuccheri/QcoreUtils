"""
Create input string for Qcore app test of 'symmetry_reduction'
for available Bravais lattices in xTB
"""

from typing import Optional

from src.utils import Set
from src import qcore_input_strings as qcore_input


def xtb_symmetry_test_string(
        crystal:      Optional[dict]=None,
        sub_commands: Optional[dict]=None,
        options:      Optional[dict]=None,
        assertions:   Optional[dict]=None,
        named_result: Optional[str]=None,
        comments='') -> str:

    #keys = [key for key in options.keys()]

    options['symmetry_reduction'] = Set(False)
    input_no_symmetry = qcore_input.xtb_input_string_cleaner(
        crystal=crystal,
        options=options,
        assertions=assertions,
        named_result=named_result + "_no_symmetry",
        comments=comments)

    options['symmetry_reduction'] = Set(True)
    input_symmetry = qcore_input.xtb_input_string_cleaner(
        crystal=crystal,
        options=options,
        assertions=assertions,
        named_result=named_result + "_no_symmetry",
        comments=comments)

    return input_no_symmetry + input_symmetry

