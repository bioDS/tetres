import os
import subprocess


def run_treeannotator(tree_file, outfile, beast_treeannotator, burnin=10, heights="ca"):
    if not os.path.exists(beast_treeannotator):
        raise ValueError("Path to beast-applauncher does not exist!")

    cmd = [beast_treeannotator,
           "-burnin", f"{burnin}",
           "-heights", heights,
           tree_file,
           outfile]
    subprocess.run(cmd, check=True)
