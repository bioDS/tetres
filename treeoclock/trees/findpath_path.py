import functools
import ete3
from ctypes import POINTER


from treeoclock.trees.time_tree import TimeTree
from treeoclock.trees._converter import ete3_to_ctree
from treeoclock.trees._ctrees import TREE, TREE_LIST
from treeoclock.trees.time_tree import lib


@functools.singledispatch
def findpath_path(arg):
    raise TypeError(type(arg) + " not supported.")


@findpath_path.register(TREE)
def _(t1, t2):
    lib.return_findpath.argtypes = [POINTER(TREE), POINTER(TREE)]
    lib.return_findpath.restype = TREE_LIST
    path = lib.return_findpath(t1, t2)
    return [path.trees[i] for i in range(path.num_trees)]


@findpath_path.register(ete3.Tree)
def _(t1, t2):
    lib.return_findpath.argtypes = [POINTER(TREE), POINTER(TREE)]
    lib.return_findpath.restype = TREE_LIST
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    path = lib.return_findpath(ct1, ct2)
    return [path.trees[i] for i in range(path.num_trees)]


@findpath_path.register(TimeTree)
def _(t1, t2):
    lib.return_findpath.argtypes = [POINTER(TREE), POINTER(TREE)]
    lib.return_findpath.restype = TREE_LIST
    path = lib.return_findpath(t1.ctree, t2.ctree)
    return [path.trees[i] for i in range(path.num_trees)]
