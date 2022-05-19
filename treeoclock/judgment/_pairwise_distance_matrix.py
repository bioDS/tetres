from ctypes import CDLL, POINTER, c_long
from treeoclock.trees._ctrees import TREE_LIST, TREE
import os
import numpy as np

lib = CDLL(f"{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/trees/findpath.so")


# todo testing required!
def calc_pw_distances(trees):
    # todo add a n_cores parameter to c function?
    n = len(trees)
    lib.pairwise_distance.argtypes = [POINTER(TREE_LIST), c_long, POINTER((c_long * n) * n)]
    distances = np.zeros((n, n), dtype=c_long)
    ctreelist = TREE_LIST(n, (TREE * n)(*[t.ctree for t in trees]))
    lib.pairwise_distance(ctreelist, n, distances.ctypes.data_as(POINTER((c_long * n) * n)))
    return distances


# todo outsource to new file!
def calc_pw_distances_two_sets(trees1, trees2):
    # todo add a n_cores parameter to c function?
    if len(trees1) != len(trees2):
        raise ValueError("Treesets given MUST have same length!")
    n = len(trees1)
    lib.pairwise_distance_two_sets.argtypes = [POINTER(TREE_LIST), POINTER(TREE_LIST), c_long, POINTER((c_long * n) * n)]
    distances = np.zeros((n, n), dtype=c_long)
    ctreelist1 = TREE_LIST(n, (TREE * n)(*[t.ctree for t in trees1]))
    ctreelist2 = TREE_LIST(n, (TREE * n)(*[t.ctree for t in trees2]))
    lib.pairwise_distance_two_sets(ctreelist1, ctreelist2, n, distances.ctypes.data_as(POINTER((c_long * n) * n)))
    return distances