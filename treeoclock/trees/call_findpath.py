__author__ = 'Lena Collienne and Jordan Kettles'

import os
from time_trees import ete3_to_ctree, TREE
from ctypes import POINTER, CDLL, Structure, c_char_p, c_int

# TODO temporary imports
import line_profiler
import sys

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')


# todo def smart_findpath(): decorator
# TODO decorator for findpath_distance to check type of given arguments and return the
#  correct function

def findpath_distance(t1, t2, c=False):
    # This being called with two ete3 trees
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]

    if c:
        return lib.findpath_distance(t1, t2)
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    return lib.findpath_distance(ct1, ct2)


class TREE_LIST(Structure):
    _fields_ = [('num_trees', c_int), ('trees', POINTER(TREE))]

    # def __init_(self, num_trees, trees):
    #     self.num_trees = num_trees
    #     self.trees = trees


# C function tree_to_string (for testing purposes)
tree_to_cluster_string = lib.tree_to_string
tree_to_cluster_string.argtypes = [POINTER(TREE)]
tree_to_cluster_string.restype = c_char_p

# C function return_findpath
findpath_path = lib.return_findpath
findpath_path.argtypes = [POINTER(TREE), POINTER(TREE)]
findpath_path.restype = TREE_LIST
