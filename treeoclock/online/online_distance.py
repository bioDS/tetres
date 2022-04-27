from treeoclock.trees.time_trees import TimeTree, TimeTreeSet
from treeoclock.trees._ctrees import TREE

from treeoclock.online._cstructures import PRECOMP, UPDATE

from ctypes import CDLL, POINTER, c_long
import os


online_lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/online_rnni.so')


def online_distance_neighbourhood(s: TimeTree, t: TimeTree):
    # Function to compute the distance of all neighbours of s to t via the online update

    online_lib.precomp.argtypes = [POINTER(TREE), POINTER(TREE)]
    online_lib.precomp.restype = POINTER(PRECOMP)

    online_lib.initialise_precomp.argtypes = [c_long]
    online_lib.initialise_precomp.restype = POINTER(PRECOMP)

    online_lib.precomp(s.ctree, t.ctree)
    return 0
