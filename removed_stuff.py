


# Made more sense to write a generic function I can use replace on
# def entos_structure_string(crystal, fractional=True, position_precison=5):
#     if fractional:
#         position_key = 'fractional'
#     else:
#         position_key = 'xyz'
#
#     positions = crystal[position_key]
#     species = crystal['species']
#     n_atoms = len(positions)
#     assert len(species) == len(positions)
#
#     # First line
#     first_position = list_to_string(positions.pop(0), precision=position_precison)
#     structure_string = "structure( " + position_key + "=([['" + species.pop(0) + "', " + first_position + "],\n"
#
#     indent = ' ' * (structure_string.find('[') + 1)
#
#     if n_atoms == 1:
#         structure_string = structure_string.rstrip('\n')[:-1]
#         print(structure_string)
#         return structure_string + "])\n"
#
#     if n_atoms > 2:
#         for ia,position in enumerate(positions[:-1]):
#             position = list_to_string(position, precision=5)
#             structure_string += indent + "['" + species[ia] + "', " + position + "],\n"
#
#     # End line
#     last_position = list_to_string(positions[-1], precision=5)
#     structure_string += indent + "['" + species[-1] + "', " + last_position + "]])\n"
#
#     return structure_string
#
# structure_string = entos_structure_string(pbs, fractional=True)

