import functools
import os
import re

import ete3
from ctypes import POINTER, CDLL, c_long, c_int

from treeoclock.trees._converter import ete3_to_ctree, ctree_to_ete3
from treeoclock.trees._ctrees import TREE, TREE_LIST

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')


class TimeTree:
    def __init__(self, nwk, f=0):
        self.etree = ete3.Tree(nwk, format=f)
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
    
    def neighbours(self):
        return neighbourhood(self)
    
    def rank_neighbours(self):
        return get_rank_neighbours(self)
    
    def nni_neighbours(self):
        return get_nni_neighbours(self)


def neighbourhood(tree: TimeTree):
    lib.rank_move.argtypes = [POINTER(TREE), c_long]
    lib.nni_move.argtypes = [POINTER(TREE), c_long, c_int]
    num_leaves = len(tree)
    neighbours = []

    for rank in range(num_leaves, num_leaves + num_leaves - 2):
        # Only do a rank move if it is possible!
        if not tree.ctree[rank].parent == rank + 1:
            working_copy = tree.copy()
            lib.rank_move(working_copy.ctree, rank)
            neighbours.append(TimeTree(ctree_to_ete3(working_copy.ctree).write(format=5)))
            if findpath_distance(neighbours[-1], tree) != 1:
                raise ValueError("Distance of a neighbour is supposed to be 1, something went wrong!")
        elif tree.ctree[rank].parent == rank + 1:
            for i in [0, 1]:
                working_copy = tree.copy()
                lib.nni_move(working_copy.ctree, rank, i)
                neighbours.append(TimeTree(ctree_to_ete3(working_copy.ctree).write(format=5)))
                if findpath_distance(neighbours[-1], tree) != 1:
                    raise ValueError("Distance of a neighbour is supposed to be 1, something went wrong!")

    return neighbours


def get_rank_neighbours(tree: TimeTree):
    lib.rank_move.argtypes = [POINTER(TREE), c_long]
    num_leaves = len(tree)
    rank_neighbours = []
    for rank in range(num_leaves, num_leaves + num_leaves - 2):
        # Only do a rank move if it is possible!
        if not tree.ctree[rank].parent == rank + 1:
            working_copy = tree.copy()
            lib.rank_move(working_copy.ctree, rank)
            rank_neighbours.append(TimeTree(ctree_to_ete3(working_copy.ctree).write(format=5)))
            if findpath_distance(rank_neighbours[-1], tree) != 1:
                raise ValueError("Distance of a neighbour is supposed to be 1, something went wrong!")
    return rank_neighbours


def get_nni_neighbours(tree: TimeTree):
    lib.nni_move.argtypes = [POINTER(TREE), c_long, c_int]
    num_leaves = len(tree)
    nni_neighbours = []
    for rank in range(num_leaves, num_leaves + num_leaves - 2):
        if tree.ctree[rank].parent == rank + 1:
            for i in [0, 1]:
                working_copy = tree.copy()
                lib.nni_move(working_copy.ctree, rank, i)
                nni_neighbours.append(TimeTree(ctree_to_ete3(working_copy.ctree).write(format=5)))
                if findpath_distance(nni_neighbours[-1], tree) != 1:
                    raise ValueError("Distance of a neighbour is supposed to be 1, something went wrong!")
    return nni_neighbours


class TimeTreeSet:
    def __init__(self, file=None):
        if file is not None:
            self.map = get_mapping_dict(file)
            self.trees = read_nexus(file)
        else:
            self.map = {}
            self.trees = []

    def __getitem__(self, index: int):
        return self.trees[index]

    def __len__(self):
        return len(self.trees)

    def fp_distance(self, i: int, j: int):
        return findpath_distance(self.trees[i].ctree, self.trees[j].ctree)

    def fp_path(self, i: int, j: int):
        return findpath_path(self.trees[i].ctree, self.trees[j].ctree)


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

    tt = TimeTree("((1:3,5:3):1,(4:2,(3:1,2:1):1):2);")

    tt.etree.render(f'/Users/larsberling/Desktop/CodingMA/test.png')

    # Defining a layout to display internal node names in the plot
    def my_layout(node):
        if node.is_leaf():
            # If terminal node, draws its name
            name_face = ete3.AttrFace("name")
        else:
            # If internal node, draws label with smaller font size
            name_face = ete3.AttrFace("name", fsize=10)
        # Adds the name face to the image at the preferred position
        ete3.faces.add_face_to_node(name_face, node, column=0, position="branch-right")


    ts = ete3.TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = my_layout
    ts.show_branch_length = True
    ts.show_scale = False


    tt.etree.show(tree_style=ts)

    d_name = 'Dengue'
    # myts = TimeTree("((1:3,5:3):1,(4:2,(3:1,2:1):1):2);")

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')
    # print(myts[0].get_newick())
    # print(ctree_to_ete3(myts[0].ctree).write(format=3))

    # def my_layout(node):
    #     if node.is_leaf():
    #         # If terminal node, draws its name
    #         name_face = ete3.AttrFace("name")
    #     else:
    #         # If internal node, draws label with smaller font size
    #         name_face = ete3.AttrFace("name", fsize=10)
    #     # Adds the name face to the image at the preferred position
    #     ete3.faces.add_face_to_node(name_face, node, column=0, position="branch-right")
    #
    #
    # # n[0].render(f'/Users/larsberling/Desktop/CodingMA/test.png')
    # ts = ete3.TreeStyle()
    # ts.show_leaf_name = False
    # ts.layout_fn = my_layout
    # ts.show_branch_length = True
    # ts.show_scale = False
    # n[0].show(tree_style=ts)


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

