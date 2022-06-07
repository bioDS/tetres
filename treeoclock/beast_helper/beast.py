import os


def run_beast(xml_file, beast):
    if not os.path.exists(beast):
        raise ValueError("Path to beast does not exist!")
    # os.chdir(os.path.dirname(xml_file))
    os.system(f"{beast} -beagle_cpu -working -overwrite -threads -1 {xml_file} > {os.path.dirname(xml_file)}/command_output.txt 2>&1")
    return 0
