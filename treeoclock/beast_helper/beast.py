import os
import sys


def run_beast(xml_file, beast, chains=4):
    if not os.path.exists(beast):
        raise ValueError("Path to beast does not exist!")
    if type(chains) is not int or chains < 0:
        raise ValueError("Chains need to be an int > 0!")
    cwd = os.getcwd()  # Getting current Working directory
    for i in range(1, chains+1):
        try:
            os.mkdir(f"{os.path.dirname(xml_file)}/run{i}")
        except FileExistsError:
            sys.exit(f"Folder run{i} already exists, won't overwrite!")
        os.chdir(f"{os.path.dirname(xml_file)}/run{i}")
        os.system(f"{beast} -beagle_cpu -overwrite -threads -1 {xml_file} > {os.path.dirname(xml_file)}/run{i}/command_output.txt 2>&1")
    os.chdir(cwd)  # Resetting the working directory
    return 0
