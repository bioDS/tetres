__author__ = 'Lena Collienne'

import os
from ctypes import *

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')

class NODE(Structure):
    _fields_ = [('parent', c_long), ('children', c_long * 2)] # The order of arguments here matters! Needs to be the same as in C code!


class TREE(Structure):
    _fields_ =  [('num_leaves', c_long), ('tree', POINTER(NODE))] # Everything from struct definition in C

    def __init_(self, num_leaves, tree, id):
        self.num_leaves = num_leaves
        self.tree = tree


# C function tree_to_string for testing purposes
py_tts = lib.tree_to_string
py_tts.argtypes= [POINTER(TREE)]
py_tts.restype = c_char_p
