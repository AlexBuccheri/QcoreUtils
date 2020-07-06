""" See if Ewald converges """

import collections
import numpy as np
import subprocess
import json

from crystal_system import cubic
from pymatgen_wrappers import cif_parser_wrapper
import qcore_input_strings as qcore_input
from utils import Set


def convention_sodium_chloride(named_result, ewald):
    fname = '../' + cubic.conventional_fcc_cifs['sodium_chloride'].file
    crystal = cif_parser_wrapper(fname, is_primitive_cell=False, fractional=True, bravais='cubic')
    settings = collections.OrderedDict([
        ('h0_cutoff',               Set(40, 'bohr')),
        ('overlap_cutoff',          Set(40, 'bohr')),
        ('repulsive_cutoff',        Set(40, 'bohr')),
        ('ewald_real_cutoff',       Set(ewald['real'], 'bohr')),
        ('ewald_reciprocal_cutoff', Set(ewald['reciprocal'])),
        ('ewald_alpha',             Set(ewald['alpha'])),
        ('monkhorst_pack',          Set([4, 4, 4])),
        ('symmetry_reduction',      Set(True)),
        ('temperature',             Set(0, 'kelvin'))
    ])

    input_string = qcore_input.xtb_input_string(crystal, settings, named_result=named_result)
    return input_string

named_result = "conventional_nacl"

def run_qcore(input_string):
    qcore_exe = '/Users/alexanderbuccheri/Codes/entos/cmake-build-debug/entos'
    qcore_command = [qcore_exe, '--format', 'json', '-s', input_string.replace('\n', ' ')]
    try:
        # Can't write any stderr to stdout as this will mess up the JSON format and hence can't parse
        qcore_json_result = subprocess.check_output(qcore_command, stderr=subprocess.DEVNULL) #, stderr=subprocess.STDOUT).decode("utf-8")
        return json.loads(qcore_json_result)
    except subprocess.CalledProcessError:  # as error:
        #print("subprocess error:", error.returncode, "found:", error.output)
        return {named_result: {"energy": 0}}




def sweep_values():
    # Should also do this for primitive NaCl
    # Want to repeat this for real = 10, 15, 25, 30
    # Want to see if this will need redoing for higher k-sampling

    ewald = {'real': 20}
    k_cutoffs = [1, 2, 3, 4, 6, 8, 10]
    alphas = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1]

    for k in k_cutoffs:
        for alpha in alphas:
            ewald['reciprocal'] = k
            ewald['alpha'] = alpha
            input_string = convention_sodium_chloride(named_result, ewald)
            result = run_qcore(input_string)
            print(ewald['real'], ewald['reciprocal'], ewald['alpha'], result[named_result]['energy'])


# Procedure described here: https://www.scd.stfc.ac.uk/Pages/DL_POLY-FAQs.aspx#FAQQ5
# Plot the Coulombic energy versus alpha. If the Ewald sum is correctly converged you will see a plateau in the plot.
# Repeat this for smaller k_max values until one can't go smaller without affecting the energy
def ewald_energy_minimum():
    ewald_real = 5 #20
    k_max = 20       #1
    alpha_mid = 3.2 / ewald_real
    # Expects items for concatenation in a tuple, hence extra brackets
    alphas = np.concatenate((np.array([0.01, 0.02, 0.05, 0.1]),
                             np.linspace(0.8* alpha_mid, 1.2*alpha_mid, num=21),
                             np.array([0.5, 1, 2])))

    for alpha in alphas:
        ewald = {'real': ewald_real, 'reciprocal':k_max, 'alpha': alpha}
        input_string = convention_sodium_chloride(named_result, ewald)
        result = run_qcore(input_string)
        print(ewald['real'], ewald['reciprocal'], ewald['alpha'], result[named_result]['energy'])


ewald_energy_minimum()