from multiprocessing.pool import ThreadPool as Pool
from multiprocessing import Pool as mpPool

from treeoclock.trees._converter import ete3_to_ctree
from treeoclock.trees.time_trees import TimeTree, TimeTreeSet, findpath_distance


def compute_sos_mt(t: TimeTree, trees: TimeTreeSet, n_cores: int = None) -> int:
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
        dists = p.starmap(findpath_distance, [(t.ctree, i.ctree) for i in trees])
    return sum(i * i for i in dists)


def compute_sos(t: TimeTree, trees: TimeTreeSet):
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
        sos += findpath_distance(t.ctree, i.ctree)**2
    return sos


# def findpath_distance_mp(t1, t2):
#     return findpath_distance(ete3_to_ctree(t1), ete3_to_ctree(t2))


import ctypes, os
from ctypes import CDLL, POINTER
from treeoclock.trees._ctrees import TREE_LIST, TREE


lib = CDLL(f"{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/trees/findpath.so")


def compute_sos_omp(t: TimeTree, trees: TimeTreeSet, n_cores: int = 2) -> int:
    lib.sum_of_squares.argtypes = [POINTER(TREE), POINTER(TREE_LIST), ctypes.c_int]
    lib.sum_of_squares.restype = ctypes.c_long

    num_trees = len(trees)
    ctreelist = TREE_LIST(num_trees, (TREE * num_trees)(*[t.ctree for t in trees]))

    sos = lib.sum_of_squares(t.ctree, ctreelist, n_cores)
    return sos


if __name__ == '__main__':
    d_name = "dispg2d_australia_small_active_1"
    # d_name = "binary_single_cell_K047_gamma_beta.667"
    # d_name = "mascot_results_newzealand_trees_COVID19_MASCOT_large_time"
    # d_name = 'RSV2'
    
    from timeit import default_timer as timer
    from random import randint

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')
    print(len(myts), len(myts[0]))

    times = []
    times_mp = []
    
    for _ in range(10):
        index = randint(0, len(myts)-1)
        s = timer()
        sos1 = compute_sos_omp(myts[index], myts, n_cores=-1)
        times.append(timer() - s)

        s = timer()
        sos2 = compute_sos_mt(myts[index], myts, n_cores=None)
        times_mp.append(timer() - s)
        
        if not sos1 == sos2:
            print("FAIL", sos1, sos2)

    import numpy as np
    print(np.mean(times), np.mean(times_mp))
