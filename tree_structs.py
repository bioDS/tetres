__author__ = 'Lena Collienne and Jordan Kettles'

import os
from ctypes import *
import numpy as np
from numpy.ctypeslib import ndpointer

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')

class NODE(Structure):
    _fields_ = [('parent', c_long), ('children', c_long * 2), ('time', c_long)] # The order of arguments here matters! Needs to be the same as in C code!
    def __init_(self, parent, children, time):
        self.parent = parent
        self.children = children
        self.time = time


class TREE(Structure):
    _fields_ =  [('num_leaves', c_long), ('tree', POINTER(NODE)), ('root_time', c_long)] # Everything from struct definition in C
    def __init_(self, num_leaves, tree, root_time):
        self.num_leaves = num_leaves
        self.tree = tree
        self.root_time = root_time


class TREE_LIST(Structure):
    _fields_ = [('num_trees', c_int), ('trees', POINTER(TREE))]
    def __init_(self, num_trees, trees):
        self.num_trees = num_trees
        self.trees = trees



# C function tree_to_string (for testing purposes)
tree_to_cluster_string = lib.tree_to_string
tree_to_cluster_string.argtypes = [POINTER(TREE)]
tree_to_cluster_string.restype = c_char_p


# C function findpath_distance
findpath_distance = lib.findpath_distance
findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
findpath_distance.restype = c_long

# C function return_findpath
findpath_path = lib.return_findpath
findpath_path.argtypes = [POINTER(TREE), POINTER(TREE)]
findpath_path.restype = TREE_LIST
