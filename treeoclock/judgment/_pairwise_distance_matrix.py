from ctypes import CDLL, POINTER, c_long
from treeoclock.trees._ctrees import TREE_LIST, TREE
import os
import numpy as np
from multiprocessing import Pool
import itertools
from share_array.share_array import get_shared_array, make_shared_array


lib = CDLL(f"{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/trees/findpath.so")


def _rf(t1, t2, i, j):
    matrix = get_shared_array('distances')
    matrix[i, j] = t1.robinson_foulds(t2)[0]
    return 0


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
        make_shared_array(distances, name='distances')  # create shared memory array from numpy array
        # shared_array = get_shared_array('distances')  # get shared memory array as numpy array
        with Pool(None) as p:
            p.starmap(_rf, [(etrees[i], etrees[j], i, j) for i,j in itertools.combinations(range(n), 2)])
        distances =  get_shared_array('distances')
        try:
            del globals()['distances']
        except KeyError:
            pass
        return distances


def calc_pw_distances_two_sets(trees1, trees2, rf: bool = False):
    # todo add a n_cores parameter to c function?
    n1 = len(trees1)
    n2 = len(trees2)
    if not rf:
        lib.pairwise_distance_two_sets.argtypes = [POINTER(TREE_LIST), POINTER(TREE_LIST), c_long, c_long, POINTER((c_long * n1) * n2)]
        distances = np.zeros((n1, n2), dtype=c_long)
        ctreelist1 = TREE_LIST(n1, (TREE * n1)(*[t.ctree for t in trees1]))
        ctreelist2 = TREE_LIST(n2, (TREE * n2)(*[t.ctree for t in trees2]))
        lib.pairwise_distance_two_sets(ctreelist1, ctreelist2, n1, n2, distances.ctypes.data_as(POINTER((c_long * n1) * n2)))
        return distances
    else:
        distances = np.zeros((n1, n2), dtype=c_long)
        etrees1 = [t.etree for t in trees1]
        etrees2 = [t.etree for t in trees2]
        make_shared_array(distances, name='distances')  # create shared memory array from numpy array
        with Pool(None) as p:
            p.starmap(_rf, [(etrees1[i], etrees2[j], i, j) for i in range(n1) for j in range(n2)])
        distances =  get_shared_array('distances')
        try:
            del globals()['distances']
        except KeyError:
            pass
        return distances
