from multiprocessing import Pool

from ..trees.time_trees import TimeTree, TimeTreeSet, findpath_distance


def compute_sos(t: TimeTree, trees: TimeTreeSet, n_cores: int) -> int:
    """
    Computes the sum of squared distances for the tree t and the set trees, using n_cores processing cores

    :param t: A tree
    :type t: TimeTree
    :param trees: A set of trees
    :type trees: TimeTreeSet
    :param n_cores: Number of processing cores
    :type n_cores: int
    :return: The sum of squared distances
    :rtype: int
    """
    with Pool(n_cores) as p:
        dists = p.starmap(findpath_distance, [(t.etree, i.etree) for i in trees])
    return sum(i * i for i in dists)
