from multiprocessing import Pool

from ..trees.time_trees import TimeTree, TimeTreeSet, findpath_distance


def compute_sos(t: TimeTree, trees: TimeTreeSet, n_cores: int):
    with Pool(n_cores) as p:
        dists = p.starmap(findpath_distance, [(t.ctree, i.ctree) for i in trees])
    return sum(i*i for i in dists)

