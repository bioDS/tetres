from ctypes import CDLL, POINTER, c_long
from treeoclock.trees._ctrees import TREE_LIST, TREE
import os
import numpy as np
from multiprocessing import Pool
import itertools


lib = CDLL(f"{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/trees/findpath.so")


def _rf(t1, t2, i, j, matrix):
    matrix[i, j] = t1.robinson_foulds(t2)[0]
    return 0


# todo testing required!
def calc_pw_distances(trees, rf: bool = False):
    # todo add a n_cores parameter to c function?
    n = len(trees)
    if not rf:
        lib.pairwise_distance.argtypes = [POINTER(TREE_LIST), c_long, POINTER((c_long * n) * n)]
        distances = np.zeros((n, n), dtype=c_long)
        ctreelist = TREE_LIST(n, (TREE * n)(*[t.ctree for t in trees]))
        lib.pairwise_distance(ctreelist, n, distances.ctypes.data_as(POINTER((c_long * n) * n)))
        return distances
    else:
        distances = np.zeros((n, n), dtype=c_long)
        etrees = [t.etree for t in trees]
        with Pool(None) as p:
            p.starmap(_rf, [(etrees[i], etrees[j], i, j, distances) for i,j in itertools.combinations(range(n), 2)])
        return distances


# todo testing required
def calc_pw_distances_two_sets(trees1, trees2, rf: bool = False):
    # todo add a n_cores parameter to c function?
    if len(trees1) != len(trees2):
        raise ValueError("Treesets given MUST have same length!")
    n = len(trees1)
    if not rf:
        lib.pairwise_distance_two_sets.argtypes = [POINTER(TREE_LIST), POINTER(TREE_LIST), c_long, POINTER((c_long * n) * n)]
        distances = np.zeros((n, n), dtype=c_long)
        ctreelist1 = TREE_LIST(n, (TREE * n)(*[t.ctree for t in trees1]))
        ctreelist2 = TREE_LIST(n, (TREE * n)(*[t.ctree for t in trees2]))
        lib.pairwise_distance_two_sets(ctreelist1, ctreelist2, n, distances.ctypes.data_as(POINTER((c_long * n) * n)))
        return distances
    else:
        distances = np.zeros((n, n), dtype=c_long)
        etrees1 = [t.etree for t in trees1]
        etrees2 = [t.etree for t in trees2]

        with Pool(None) as p:
            p.starmap(_rf, [(etrees1[i], etrees2[j], i, j, distances) for i, j in
                            itertools.product(range(n), repeat=2)])
        return distances