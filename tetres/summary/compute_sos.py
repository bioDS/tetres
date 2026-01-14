from multiprocessing.pool import ThreadPool as Pool
import ctypes, os
from ctypes import CDLL, POINTER

from tetres.trees._ctrees import TREE_LIST, TREE
from tetres.trees.time_trees import TimeTree, TimeTreeSet, findpath_distance


def compute_sos_mt(t: TimeTree, trees: TimeTreeSet, n_cores: int = None, norm=False) -> int:
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
        dists = p.starmap(findpath_distance, [(t.ctree, i.ctree, norm) for i in trees])
    return sum(i * i for i in dists)


def compute_sos(t: TimeTree, trees: TimeTreeSet, norm=False):
    """
    Computes the sum of squared distances for the tree t and the set trees

    :param t: A tree
    :type t: TimeTree
    :param trees: A set of trees
    :type trees: TimeTreeSet
    :return: The sum of squared distances
    :rtype: int
    """
    sos = 0
    for i in trees:
        sos += findpath_distance(t.ctree, i.ctree, norm=norm)**2
    return sos


def compute_sos_omp(t: TimeTree, trees: TimeTreeSet, n_cores: int = None) -> int:
    lib = CDLL(f"{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/trees/findpath.so")
    lib.sum_of_squares.argtypes = [POINTER(TREE), POINTER(TREE_LIST), ctypes.c_int]
    lib.sum_of_squares.restype = ctypes.c_long

    num_trees = len(trees)
    ctreelist = TREE_LIST(num_trees, (TREE * num_trees)(*[t.ctree for t in trees]))

    if n_cores is None:
        # c does not allow None, -1 is the argument for all possible cores
        n_cores = -1

    sos = lib.sum_of_squares(t.ctree, ctreelist, n_cores)
    return sos


from multiprocessing import get_context
# Global Pool
_POOL = None
_T = None
_NB = None


def _init_pool(t, nb_taxa):
    """Initializer for Pool: set global references."""
    global _T, _NB
    _T = t
    _NB = nb_taxa


def _worker_chunk(trees_chunk):
    """Top-level function for Pool mapping."""
    from GelmanRubin import __hop_min_distance
    results = [__hop_min_distance(_T, tree, _NB) for tree in trees_chunk]
    return results


def compute_hop_sos_mt(t, trees, n_cores: int = None) -> int:
    from GelmanRubin import __hop_min_distance
    from TreeVec import get_nb_taxa
    import numpy as np

    global _POOL, _T, _NB

    if n_cores is None:
        import os
        n_cores = os.cpu_count()

    _T = t
    _NB = get_nb_taxa(t[0])

    if _POOL is None:
        _POOL = get_context("fork").Pool(processes=n_cores, initializer=_init_pool,
                                         initargs=(t, _NB))

    # Split trees into chunks
    chunks = np.array_split(trees, n_cores)

    # Map each chunk to a worker
    dists_chunks = _POOL.map(_worker_chunk, chunks)

    # Flatten list of lists
    dists = [dist for chunk in dists_chunks for dist in chunk]

    return np.sum(np.square(dists))
