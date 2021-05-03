__author__ = 'Lena Collienne and Jordan Kettles'

import os
from time_trees import ete3_to_ctree, TREE
from ctypes import POINTER, CDLL, Structure, c_char_p, c_int


class TREE_LIST(Structure):
    _fields_ = [('num_trees', c_int), ('trees', POINTER(TREE))]

    def __init_(self, num_trees, trees):
        self.num_trees = num_trees
        self.trees = trees


# C function tree_to_string (for testing purposes)
tree_to_cluster_string = lib.tree_to_string
tree_to_cluster_string.argtypes = [POINTER(TREE)]
tree_to_cluster_string.restype = c_char_p

# C function return_findpath
findpath_path = lib.return_findpath
findpath_path.argtypes = [POINTER(TREE), POINTER(TREE)]
findpath_path.restype = TREE_LIST
