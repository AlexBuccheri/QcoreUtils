""" Module to run qcore"""

import subprocess
import json

# Specific to Alex's machine
_PATH_TO_ENTOS = "/Users/alexanderbuccheri/Codes/entos/"
_QCORE_EXES = {"debug": _PATH_TO_ENTOS + 'cmake-build-debug/qcore'}


def get_named_result(input_string: str) -> str:
    """ Get first named result in qcore input string """
    index = input_string.find(":")
    if index != -1:
        return input_string[:index]
    else:
        raise Exception("Unable to find named result in ", input_string)


#TODO(Alex) Tidy up the return if error
def run_qcore(input_string: str, exe_type="debug") -> dict:
    """
    Run qcore, get the JSON output via stdout and
    return the results dictionary

    """
    qcore_exe = _QCORE_EXES[exe_type]
    named_result = get_named_result(input_string)
    qcore_command = [qcore_exe, '--format', 'json', '-s', input_string.replace('\n', ' ')]
    try:
        # Can't write any stderr to stdout as this will mess up the JSON format and hence can't parse
        qcore_json_result = subprocess.check_output(qcore_command, stderr=subprocess.DEVNULL) #, stderr=subprocess.STDOUT).decode("utf-8")
        return json.loads(qcore_json_result)
    except subprocess.CalledProcessError:  # as error:
        #print("subprocess error:", error.returncode, "found:", error.output)
        return {named_result: {}}