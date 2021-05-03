__author__ = 'Lena Collienne, Lars Berling, Jordan Kettles'

import re
import os
import ete3
from ctypes import c_long, Structure, POINTER, CDLL, c_int

# TODO temporary imports
import line_profiler

# TODO doumentation for the classes missing

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')


class NODE(Structure):
    _fields_ = [('parent', c_long), ('children', c_long * 2),
                ('time', c_long)]  # The order of arguments here matters! Needs to be the same as in C code!

    def __init_(self, parent, children, time):
        self.parent = -1
        self.children = [-1, -1]
        self.time = 0


class TREE(Structure):
    _fields_ = [('num_leaves', c_long), ('tree', POINTER(NODE)),
                ('root_time', c_long)]  # Everything from struct definition in C

    def __init_(self, num_leaves, tree, root_time):
        self.num_leaves = num_leaves
        self.tree = tree
        self.root_time = root_time


class TREE_LIST(Structure):
    _fields_ = [('num_trees', c_int), ('trees', POINTER(TREE))]

    def __init_(self, num_trees, trees):
        self.num_trees = num_trees
        self.trees = trees


class TimeTree:
    def __init__(self, nwk):
        self.etree = ete3.Tree(nwk)
        self.ctree = ete3_to_ctree(self.etree)

    def fp_distance(self, tree):
        return findpath_distance(self.ctree, tree.ctree, True)

    def get_newick(self):
        return self.etree.write(format=5)

    def write_newick(self, file):
        return self.etree.write(file, format=5)

    # TODO one_neighbourhood function



def findpath_path(t1, t2, c=False):
    # C function return_findpath
    lib.return_findpath.argtypes = [POINTER(TREE), POINTER(TREE)]
    lib.return_findpath.restype = TREE_LIST
    if c:
        return lib.return_findpath(t1, t2)
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    return lib.return_findpath(ct1, ct2)


def findpath_distance(t1, t2, c=False):
    # This being called with two ete3 trees
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]

    if c:
        return lib.findpath_distance(t1, t2)
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    return lib.findpath_distance(ct1, ct2)


def ctree_to_ete3(ctree):
    nl = ctree.num_leaves
    nn = (nl * 2) - 2  # number of nodes - 1, max index in ctree.tree

    def traverse(node):
        nonlocal ctree
        nonlocal nl

        # Curent node is an internal node
        cur_t = ete3.Tree()
        if node.children[0] >= nl:
            cur_t.add_child(traverse(ctree.tree[node.children[0]]))
        else:
            cur_t.add_child(name=node.children[0] + 1)
        if node.children[1] >= nl:
            cur_t.add_child(traverse(ctree.tree[node.children[1]]))
        else:
            cur_t.add_child(name=node.children[1] + 1)
        return cur_t

    t = traverse(ctree.tree[nn])
    return t


# @profile
def ete3_to_ctree(tree):
    distances = []
    node2leaves = tree.get_cached_content()
    # tree_root = tree.get_tree_root()

    index = 0
    for node in tree.traverse('levelorder'):
        if len(node2leaves[node]) != 1:
            # if node.children:
            # if not node.is_leaf():
            if index == 0:
                node.name = node.dist
            else:
                # node.name = node.get_distance(tree_root))  # get_distance() is a very slow function
                node.name = node.dist + node.up.name  # Top down for loop makes this possible
            index += 1
            distances.append(node.name)

    num_nodes = len(node2leaves)  # Number of nodes
    num_leaves = int(((num_nodes - 1) / 2) + 1)
    node_list = (NODE * num_nodes)()

    distances = sorted(distances)
    # if not len(distances.keys()) == num_leaves - 1:
    if not len(set(distances)) == num_leaves - 1:
        sys.exit('Distances to root not unique! \n'
                 'This has to be resolved!')
    for node in tree.traverse('levelorder'):
        if len(node2leaves[node]) == 1:
            # if not node.children:
            # if node.is_leaf():
            node_list[int(node.name) - 1].parent = num_nodes - (distances.index(node.up.name) + 1)
        else:
            if node.name == 0.0:
                node_list[num_nodes - (distances.index(node.name) + 1)].parent = num_nodes - 1
            else:
                node_list[num_nodes - (distances.index(node.name) + 1)].parent = \
                    num_nodes - (distances.index(node.up.name) + 1)
            # current_children = node.get_children()  # get_children() is slow
            current_children = node.children
            if len(current_children) != 2:
                sys.exit('Not a binary tree!')

            # Child 0
            if len(node2leaves[current_children[0]]) == 1:
                # if not current_children[0].children:
                # if current_children[0].is_leaf():
                node_list[num_nodes - (distances.index(node.name) + 1)].children[0] = \
                    int(current_children[0].name) - 1
            else:
                node_list[num_nodes - (distances.index(node.name) + 1)].children[0] = \
                    num_nodes - (distances.index(current_children[0].name) + 1)
            # Child 1
            if len(node2leaves[current_children[1]]) == 1:
                # if not current_children[1].children:
                # if current_children[1].is_leaf():
                node_list[num_nodes - (distances.index(node.name) + 1)].children[1] = \
                    int(current_children[1].name) - 1
            else:
                node_list[num_nodes - (distances.index(node.name) + 1)].children[1] = \
                    num_nodes - (distances.index(current_children[1].name) + 1)
    return TREE(num_leaves, node_list, -1)


def get_mapping_dict(file: str) -> dict:
    """
    Returns the taxon mapping of the nexus file as a dictionary

    :param file: A nexus file path
    :type file: str
    :return: Dictionary containing the mapping of taxa(values) to int(keys)
    :rtype: dict {int --> str}
    """

    begin_map = re.compile('\tTranslate\n', re.I)
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


def read_nexus(file, c=False):
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
                if c:
                    trees.append(ete3_to_ctree(ete3.Tree(re.sub(brackets, "", tree_string))))
                else:
                    trees.append(ete3.Tree(re.sub(brackets, "", tree_string)))
    # if c:
    #     return [ete3_to_ctree(tr) for tr in trees]
    return trees


def my_trees_read(file):
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


if __name__ == '__main__':

    # TODO need one_neighbourhood
    # TODO TREELIST to ete3 trees
    # TODO TREELIST extract one tree with an index ?

    import sys
    import random
    from timeit import default_timer as timer

    d_name = 'Dengue'

    t = read_nexus(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees', c=False)
    
    ct = read_nexus(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees', c=True)
    
    # myt = my_trees_read(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')
    
    p = findpath_path(ct[0], ct[1], True)
    ep = findpath_path(t[0], t[1], False)

    # from pympler import asizeof
    #
    # print(asizeof.asizeof(t))
    # print(asizeof.asizeof(ct))
    # print(asizeof.asizeof(myt))
