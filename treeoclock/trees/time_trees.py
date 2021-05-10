import functools
import os
import re

import ete3
from ctypes import POINTER, CDLL, c_long, c_int, memmove, sizeof, byref

from treeoclock.trees._converter import ete3_to_ctree, ctree_to_ete3
from treeoclock.trees._ctrees import TREE, TREE_LIST

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')


class TimeTree:
    def __init__(self, nwk):
        self.etree = ete3.Tree(nwk)
        self.ctree = ete3_to_ctree(self.etree)

    def __len__(self):
        return self.ctree.num_leaves

    def fp_distance(self, tree):
        return findpath_distance(self.ctree, tree.ctree)
    
    def fp_path(self, tree):
        return findpath_path(self.ctree, tree.ctree)

    def get_newick(self, f=5):
        return self.etree.write(format=f)
    
    def copy(self):
        return TimeTree(self.get_newick())

    # TODO one_neighbourhood function


def neighbourhood(tree: TimeTree):
    lib.rank_move.argtypes = [POINTER(TREE), c_long]

    num_leaves = len(tree)

    # All possible rank moves
    print(tree.get_newick(f=9))
    copy = tree.copy()
    for rank in range(num_leaves - 1):
        lib.rank_move(copy.ctree, rank)  # TODO Changes the tree inline, but that should not be the case!
        print(ctree_to_ete3(copy.ctree).write(format=9) == ctree_to_ete3(tree.ctree).write(format=9))

    print("------")
    print(ctree_to_ete3(tree.ctree).write(format=9))
    print(ctree_to_ete3(copy.ctree).write(format=9))
    # All NNI moves
    lib.nni_move.argtypes = [POINTER(TREE), c_long, c_int]
    # for rank in range(num_leaves - 1):
        # First logic part: does the parent have rank+1:
            # lib.nni_move

    # TODO how to check if nni move is possible for a given index i and i+1 ?
    #  Maybe get exceptions from the c code ? or do the same logic here so that the c code does not fail
    # This should return a numpy array with TimeTree objects
    return 0


class TimeTreeSet:
    def __init__(self, file):
        self.map = get_mapping_dict(file)
        self.trees = read_nexus(file)

    def __getitem__(self, index):
        return self.trees[index]

    def __len__(self):
        return len(self.trees)

    def fp_distance(self, i, j):
        return findpath_distance(self.trees[i].ctree, self.trees[j].ctree)

    def fp_path(self, i, j):
        return findpath_path(self.trees[i].ctree, self.trees[j].ctree)
    # TODO make it possible to iterate over a TimeTreeSet for i in TTS: i is then the tree


def read_nexus(file: str) -> list:
    # re_tree returns nwk string without the root height and no ; in the end
    re_tree = re.compile("\t?tree .*=? (.*$)", flags=re.I | re.MULTILINE)
    # Used to delete the ; and a potential branchlength of the root
    # name_dict = get_mapping_dict(file)  # Save tree label names in dict
    brackets = re.compile(r'\[[^\]]*\]')  # Used to delete info in []

    trees = []
    with open(file, 'r') as f:
        for line in f:
            if re_tree.match(line):
                tree_string = f'{re.split(re_tree, line)[1][:re.split(re_tree, line)[1].rfind(")") + 1]};'
                trees.append(TimeTree(re.sub(brackets, "", tree_string)))
    return trees


def get_mapping_dict(file: str) -> dict:
    """
    Returns the taxon mapping of the nexus file as a dictionary

    :param file: A nexus file path
    :type file: str
    :return: Dictionary containing the mapping of taxa(values) to int(keys)
    :rtype: dict {int --> str}
    """

    begin_map = re.compile('\t?translate\n', re.I)
    end = re.compile('\t?;\n?')

    mapping = {}

    begin = False
    with open(file) as f:
        for line in f:
            if begin:
                if end.match(line):
                    break
                split = line.split()

                mapping[int(split[0])] = split[1][:-1] if split[1][-1] == "," else split[1]

            if begin_map.match(line):
                begin = True
    return mapping


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


if __name__ == '__main__':

    d_name = 'Dengue'

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    n = neighbourhood(myts[0])

    # print(len(myts[0]))
    # print(len(myts))
    #
    # print(myts[0].fp_distance(myts[0]))
    #
    # for t in myts:
    #     print(t.get_newick())
    #
    # print(myts.map)

    # from pympler import asizeof
    #
    # print(asizeof.asizeof(t))
    # print(asizeof.asizeof(ct))
    # print(asizeof.asizeof(myt))

