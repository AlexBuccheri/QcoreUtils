"""
Resolve xtb potential
"""

import enum

def xtb_potential_str(potential_type: enum.Enum) -> str:
    """
    Set xTB potential
    :param potential_type:
    :return:
    """
    if potential_type == PotentialType.TRUNCATED:
        return 'truncated'

    elif potential_type == PotentialType.Full:
        return 'full'
    else:
        raise Exception("Invalid choice for potential_type: ", potential_type)


class PotentialType(enum.Enum):
    TRUNCATED = 1
    FULL = 2