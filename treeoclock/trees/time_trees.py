import functools
import os
import re
import warnings
import ete3

from ctypes import POINTER, CDLL, c_long, c_int
from treeoclock.trees._converter import ete3_to_ctree, ctree_to_ete3
from treeoclock.trees._ctrees import TREE, TREE_LIST

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')


class DifferentNbrTaxa(Exception):
    """Raised when two trees have differnt number of taxa (findpath call)"""
    pass


class TimeTree:
    def __init__(self, nwk, f=0):
        self.etree = ete3.Tree(nwk, format=f)
        self.ctree = ete3_to_ctree(self.etree)

    def __len__(self):
        return self.ctree.num_leaves

    def fp_distance(self, tree, norm=False):
        return findpath_distance(self.ctree, tree.ctree, norm=norm)

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

    def get_clades(self):
        return nwk_to_cluster(self.get_newick(f=9))

    def apply_new_taxa_map(self, new_map, old_map):

        if not sorted(new_map.keys()) == sorted(old_map.keys()):
            raise ValueError("New map does not fit, taxa should be encoded by integers 1,...,n!")
        if not sorted(new_map.values()) == sorted(old_map.values()):
            raise ValueError("New map does not fit, different taxa names!")
        cur_nwk = self.get_newick()
        re_taxa = re.compile('([0-9]+)([\\[:])')
        new_map_reversed = {v: k for (k, v) in new_map.items()}
        new_newick = re_taxa.sub(lambda m: m.group().replace(m.group(1),
                                                             str(new_map_reversed[
                                                                     old_map[int(m.group(1))]]))
                                 , cur_nwk)

        return TimeTree(new_newick)


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
            if os.path.exists(file):
                self.map = get_mapping_dict(file)
                self.trees = read_nexus(file)
            else:
                raise FileNotFoundError(file)
        else:
            self.map = {}
            self.trees = []
        self.common_clades = set()

    def __getitem__(self, index):
        if isinstance(index, int):
            if index >= len(self.trees):
                raise IndexError("Index out of range!")
            return self.trees[index]
        elif isinstance(index, slice):
            ret = TimeTreeSet()
            ret.map = self.map
            ret.trees = self.trees[index]
            return ret
        else:
            raise TypeError("Given Type for getitem not implemented!")

    def __len__(self):
        return len(self.trees)

    def fp_distance(self, i: int, j: int, norm=False):
        return findpath_distance(self.trees[i].ctree, self.trees[j].ctree, norm=norm)

    def fp_path(self, i: int, j: int):
        return findpath_path(self.trees[i].ctree, self.trees[j].ctree)

    def copy(self):
        return self.trees.copy()

    def get_common_clades(self):
        if len(self.common_clades) == 0:
            self.common_clades = self.trees[0].get_clades()
            for i in range(1, len(self.trees)):
                self.common_clades = set.intersection(self.common_clades, self[i].get_clades())
        return self.common_clades

    def change_mapping(self, new_map):
        if new_map == self.map:
            raise ValueError("New map is identical to old one!")
        if not sorted(self.map.keys()) == sorted(new_map.keys()):
            raise ValueError("New map does not fit, different encodings of the taxa (Should generally be 1,...,n)!")
        if not sorted(self.map.values()) == sorted(new_map.values()):
            raise ValueError("New map does not fit, different taxa names!")

        for index, t in enumerate(self.trees):
            # apply the mapping to each tree in the set!
            self.trees[index] = t.apply_new_taxa_map(new_map, self.map)  # ?
        # change the self.map in the end!
        self.map = new_map


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
def _(t1, t2, norm=False):
    if not t1.num_leaves == t2.num_leaves:
        raise DifferentNbrTaxa
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    if norm:
        return lib.findpath_distance(t1, t2)/_diam_treespace(t1.num_leaves)
    return lib.findpath_distance(t1, t2)


@findpath_distance.register(ete3.Tree)
def _(t1, t2, norm=False):
    if not len(t1.get_tree_root()) == len(t2.get_tree_root()):
        raise DifferentNbrTaxa
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    if norm:
        return lib.findpath_distance(ct1, ct2)/_diam_treespace(ct1.num_leaves)
    return lib.findpath_distance(ct1, ct2)


@findpath_distance.register(TimeTree)
def _(t1, t2, norm=False):
    if not len(t1) == len(t2):
        raise DifferentNbrTaxa
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    if norm:
        return lib.findpath_distance(t1.ctree, t2.ctree)/_diam_treespace(len(t1))
    return lib.findpath_distance(t1.ctree, t2.ctree)


@functools.singledispatch
def findpath_path(arg):
    # Function does not contain the last tree, so findpath_path(x,y)=list of trees [x, z1, z2, ...]
    raise TypeError(type(arg) + " not supported.")


# todo raising the exception if different number of taxa, different types for t1, t2


@findpath_path.register(TREE)
def _(t1, t2):
    lib.return_findpath.argtypes = [POINTER(TREE), POINTER(TREE)]
    lib.return_findpath.restype = TREE_LIST
    warnings.warn(UserWarning('The returned Tree_list object has allocated memory that needs to be freed, use lib.free_tree_list or trees.free_tree_list!'))
    return lib.return_findpath(t1, t2)


@findpath_path.register(ete3.Tree)
def _(t1, t2):
    lib.return_findpath.argtypes = [POINTER(TREE), POINTER(TREE)]
    lib.return_findpath.restype = TREE_LIST
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    warnings.warn(UserWarning('The returned Tree_list object has allocated memory that needs to be freed, use lib.free_tree_list or trees.free_tree_list!'))
    return lib.return_findpath(ct1, ct2)


@findpath_path.register(TimeTree)
def _(t1, t2):
    lib.return_findpath.argtypes = [POINTER(TREE), POINTER(TREE)]
    lib.return_findpath.restype = TREE_LIST
    warnings.warn(UserWarning('The returned Tree_list object has allocated memory that needs to be freed, use lib.free_tree_list or trees.free_tree_list!'))
    return lib.return_findpath(t1.ctree, t2.ctree)


def free_tree_list(tree_list):
    lib.free_treelist.argtypes = [TREE_LIST]
    lib.free_treelist(tree_list)
    return 0


def nwk_to_cluster(treestr):
    clades = set()
    if not (treestr[0] == '(' and treestr[-2] == ')' and treestr[-1] == ';'):
        raise Exception("Invalid tree string given! (no ';' at the end)")
    opend = []
    re_brackets = re.compile(r"\(|\)")
    for i in range(1, len(treestr) - 2):
        if treestr[i] == '(':
            opend.append(i)
        elif treestr[i] == ')':
            if not opend:
                raise Exception("Invalid tree string given! (to many ')')")
            cur = treestr[opend[-1]:i]
            clades.add(frozenset(re.sub(re_brackets, '', cur).split(',')))
            del opend[-1]
    if opend:
        raise Exception("Invalid tree string given! (to many '(')")
    return clades


def _diam_treespace(n):
    return ((n - 1) * (n - 2)) / 2
