from collections import OrderedDict
import copy
import numpy as np

from src import qcore_input_strings as qcore_input, utils
import translation_invariance.string_generator as trans

def comments_list(rotations) -> list:
    """
    Create comments list containing comment strings for each choice of rotation
    :param rotation: n x 3 array
    :return: comments list
    """

    comments = []

    for rot_patt_i in rotations:
        if np.equal(rot_patt_i, np.array([0,0,0])):
            comments.append('! No rotation')
        else:
            comments.append('! Rotation around axes: z ' + str(rot_patt_i[0]) +
                                                   ' x ' + str(rot_patt_i[1]) +
                                                   ' y ' + str(rot_patt_i[2]))
    return comments

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
        assert len(r) == 3, "Rotation must be a n x 3 matrix"
        updated_crystal = utils.rotate_lattice_vector(copy.deepcopy(crystal), r[0], r[1], r[2])
        input += trans.xtb_input_string(updated_crystal, options, assertions,
                                        named_result + '_rotation' + str(i), comments=comments[i]) + '\n'

    return input

