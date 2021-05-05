import functools
import ete3
from ctypes import POINTER


from treeoclock.trees._ctrees import TREE
from treeoclock.trees.time_trees import lib
from treeoclock.trees._converter import ete3_to_ctree
from treeoclock.trees.time_trees import TimeTree


@functools.singledispatch
def findpath_distance(arg):
    raise TypeError(type(arg) + " not supported.")


@findpath_distance.register(TREE)
def _(t1, t2):
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    return lib.findpath_distance(t1, t2)


@findpath_distance.register(ete3.Tree)
def _(t1, t2):
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    return lib.findpath_distance(ct1, ct2)


@findpath_distance.register(TimeTree)
def _(t1, t2):
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    return lib.findpath_distance(t1.ctree, t2.ctree)