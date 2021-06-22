from treeoclock.trees.time_trees import TimeTreeSet, TimeTree
from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary._tree_proposals import search_neighbourhood_greedy, NoBetterNeighbourFound

import random


def greedy(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    sos = compute_sos_mt(start, trees, n_cores=n_cores)
    centroid = start

    while True:
        try:
            centroid, sos = search_neighbourhood_greedy(t=centroid, trees=trees, t_value=sos,
                                                        n_cores=n_cores, select=select)
        except NoBetterNeighbourFound:
            # This is thrown when no neighbour has a better SoS value, i.e. the loop can be stopped
            break
    return centroid, sos


def inc_sub(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    # Initializing all parameters
    sample = TimeTreeSet()
    cen = start
    sos = compute_sos_mt(start, trees, n_cores=n_cores)
    trees_copy = trees.copy()
    cur_length = len(trees_copy)

    while cur_length > 0:
        sample.trees.extend([trees_copy.pop(random.randrange(len(trees_copy)))
                             for _ in range(min(kwargs["subsample_size"], cur_length))])
        cur_length -= kwargs["subsample_size"]

        new_cen, _ = greedy(trees=sample, n_cores=n_cores, select=select, start=start)
        new_sos = compute_sos_mt(cen, trees, n_cores=n_cores)

        if new_sos <= sos:
            sos = new_sos
            cen = new_cen
        else:
            break
    return cen, sos


def iter_sub(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    # Initializing all parameters
    sample = TimeTreeSet()
    cen = start
    sos = compute_sos_mt(start, trees, n_cores=n_cores)
    tree_len = len(trees)

    while True:
        sample.trees = [trees[(random.randrange(len(trees)))]
                        for _ in range(min(kwargs["subsample_size"], tree_len))]

        new_cen, _ = greedy(trees=sample, n_cores=n_cores, select=select, start=start)
        new_sos = compute_sos_mt(cen, trees, n_cores=n_cores)

        if new_sos <= sos:
            sos = new_sos
            cen = new_cen
        else:
            break
    return cen, sos
