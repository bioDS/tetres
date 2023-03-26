import linecache
import re


def _extract_cutoff(multichain, i, start, end, tree_file, log_file):

    if multichain.log_files[i] is not None:
        with open(log_file, "x") as f:
            f.write("\t".join(v for v in list(multichain[i].log_data.keys())))
            f.write("\n")

    with open(tree_file, "x") as f:
        f.write(f"#NEXUS\n\nBegin taxa;\n\tDimensions ntax={multichain[i].trees[i].ctree.num_leaves};\n\t\tTaxlabels\n")
        for taxa in range(1, multichain[i].trees[0].ctree.num_leaves + 1):
            f.write(f"\t\t\t{multichain[i].trees.map[taxa]}\n")
        f.write("\t\t\t;\nEnd;\nBegin trees;\n\tTranslate\n")
        for taxa in range(1, multichain[i].trees[0].ctree.num_leaves):
            f.write(f"\t\t\t{taxa} {multichain[i].trees.map[taxa]},\n")
        f.write(
            f"\t\t\t{multichain[i].trees[0].ctree.num_leaves} {multichain[i].trees.map[multichain[i].trees[0].ctree.num_leaves]}\n")
        f.write(";\n")

    sample = 1
    chain_treefile = f"{multichain.working_dir}/{multichain.tree_files[i]}"
    re_tree = re.compile("\t?tree .*=? (.*$)", flags=re.I | re.MULTILINE)
    for index in range(start, end + 1):
        offset = 10 + (2 * multichain[i].trees[0].ctree.num_leaves)

        line = index + 1 + offset
        cur_tree = linecache.getline(chain_treefile, line)
        cur_tree = f'{re.split(re_tree, cur_tree)[1][:re.split(re_tree, cur_tree)[1].rfind(")") + 1]};'
        with open(tree_file, "a") as tf:
            tf.write(f"tree STATE_{sample} = {cur_tree}\n")

        if multichain.log_files[i] is not None:
            cur_values = [v for v in multichain[i].log_data.iloc[index]]
            cur_values[0] = sample
            with open(log_file, "a") as lf:
                lf.write("\t".join([str(v) for v in cur_values]))
                lf.write("\n")
        sample += 1
    with open(tree_file, "a") as tf:
        tf.write("End;")
