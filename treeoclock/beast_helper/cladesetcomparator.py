import os


def run_cladesetcomparator(tree_file1, tree_file2, out_file_png, beast_applauncher, burnin=10):
    if not os.path.exists(beast_applauncher):
        raise ValueError("Path to beast-applauncher does not exist!")

    os.system(
        f"{beast_applauncher} CladeSetComparator "
        f"-tree1 {tree_file1} "
        f"-tree2 {tree_file2} "
        f"-png {out_file_png} "
        f"-burnin {burnin} "
        f"> {out_file_png}.txt")
