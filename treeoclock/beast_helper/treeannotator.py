import os


def run_treeannotator(tree_file, outfile, beast_treeannotator, burnin=10):
    if not os.path.exists(beast_treeannotator):
        raise ValueError("Path to beast-applauncher does not exist!")

    os.system(
        f"{beast_treeannotator} "
        f"-burnin {burnin} "
        f"-heights 'ca' "
        f"{tree_file} "
        f"{outfile}")
