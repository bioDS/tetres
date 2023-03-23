from tetres.trees.time_trees import TimeTree, TimeTreeSet
from tetres.trees._ctrees import TREE

from tetres.online._cstructures import PRECOMP, UPDATE

from ctypes import CDLL, POINTER, c_long, cast
import os, sys
import numpy as np

online_lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/online_rnni.so')


def online_distance_neighbourhood(s: TimeTree, t: TimeTree):
    # Function to compute the distance of all neighbours of s to t via the online update

    online_lib.precomp.argtypes = [POINTER(TREE), POINTER(TREE)]
    online_lib.precomp.restype = POINTER(PRECOMP)

    # online_lib.initialise_precomp.argtypes = [c_long]
    # online_lib.initialise_precomp.restype = POINTER(PRECOMP)

    online_lib.find_distance_diff.argtypes = [POINTER(TREE), POINTER(TREE), POINTER(c_long), POINTER(PRECOMP)]
    online_lib.find_distance_diff.restype = UPDATE

    precomp_data = online_lib.precomp(s.ctree, t.ctree)

    if not precomp_data.contents.dist == s.fp_distance(t):
        sys.exit("Wrong distance with precomp function!")

    num_leaves = len(s)

    # Segmentation Fault ...
    for rank in range(num_leaves, num_leaves + num_leaves - 2):
        if not s.ctree[rank].parent == rank + 1:
            online_lib.find_distance_diff(s.ctree, t.ctree, (c_long*2)(rank, 0), precomp_data)
        elif s.ctree[rank].parent == rank + 1:
            for i in [1, 2]:
                online_lib.find_distance_diff(s.ctree, t.ctree, (c_long*2)(rank, i), precomp_data)
    return 0
