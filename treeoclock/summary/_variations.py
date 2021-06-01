from treeoclock.trees.time_trees import TimeTreeSet, TimeTree

# TODO more parameters, n_cores missing
def greedy(trees: TimeTreeSet, n_cores: int, select: str):
    return len(trees)