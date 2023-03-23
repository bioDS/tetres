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
