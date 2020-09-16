# Notes
# Simple cubic: alpha-N2 with and without symmetry disagrees by ~ 1.e-5 Ha. [4,4,4] grid
# BCC:          nickel_triphosphide won't converge, Potassium 'seems ok
# but warning: solve(): system seems singular (rcond: 7.23826e-17); attempting approx solution
# so I do not know if the value is reliable
# FCC:          silicon is ok
#
# May be able to assert the total number of k-points used
# in the irreducible grid


from symmetry.crystals import *

def produce_test_string(crystal_inputs:list, header=''):
    total_input = header
    for input in crystal_inputs:
        total_input += input() + '\n'
    return total_input


header = """!Symmetry application test for all lattices for which the 
!symmetry_reduction option is available (cubic)  
!Total energy should be unaffected by utilising symmetry_reduction 
!on the k-point grid"""

all_tests = [alpha_N2, potassium, silicon, boron_nitride_hex]

print(produce_test_string(all_tests, header=header))



