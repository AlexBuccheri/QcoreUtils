from collections import OrderedDict
import copy
import numpy as np

from src import qcore_input_strings as qcore_input, utils
import translation_invariance as trans


def xtb_rotational_invariance_string(
        crystal: dict,
        options: dict,
        assertions: dict,
        rotation,
        named_result: str,
        comments=None) -> str:
    """
    Rotational invariance of the total energy means that the assertions can be the same
    """

    if comments:
        assert len(comments) == len(rotation), "len(comments) != len(rotation)"
    else:
        comments = [''] * len(rotation)

    input = ''
    for i,r in enumerate(rotation):
        assert isinstance(r, np.array), "Rotation must be a n x 3 matrix"
        assert len(r) == 3, "Rotation must be a n x 3 matrix"
        updated_crystal = utils.rotate_lattice_vector(copy.deepcopy(crystal), r)
        input += trans.xtb_input_string(updated_crystal, options, assertions,
                                        named_result + '_rotation' + str(i), comments=comments[i]) + '\n'

    return input

