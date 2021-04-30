__author__ = 'Lena Collienne and Jordan Kettles'

import os
import sys
from ctypes import *  # TODO
import line_profiler

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


# @profile
def ete3_to_ctree(tree):
    num_leaves = len(tree.get_tree_root())
    num_nodes = (num_leaves * 2) - 1
    node_list = (NODE * num_nodes)()

    # leaves = [str(i+1) for i in range(num_leaves)]

    distances = {}

    # tree_root = tree.get_tree_root()

    index = 0
    for node in tree.traverse('levelorder'):
        if node.children:
            # if not node.is_leaf():
            if index == 0:
                node.name = node.dist
            else:
                # distances.append(node.get_distance(tree_root))  # get_distance() is a very slow function
                node.name = node.dist + node.up.name  # Top down foor loop makes this possible
            distances[node.name] = index
            index += 1

    # distances = sorted(distances)
    if not len(distances.keys()) == num_leaves - 1:
        sys.exit('Distances to root not unique! \n'
                 'This has to be resolved!')
    for node in tree.traverse('levelorder'):
        if not node.children:
            # if node.is_leaf():
            node_list[int(node.name) - 1].parent = num_nodes - (distances[node.up.name] + 1)
        else:
            if node.name == 0.0:
                node_list[num_nodes - (distances[node.name] + 1)].parent = num_nodes - 1
            else:
                node_list[num_nodes - (distances[node.name] + 1)].parent = \
                    num_nodes - (distances[node.up.name] + 1)
            # current_children = node.get_children()
            current_children = node.children
            if len(current_children) != 2:
                sys.exit('Not a binary tree!')

            # Child 0
            if not current_children[0].children:
                # if current_children[0].is_leaf():
                node_list[num_nodes - (distances[node.name] + 1)].children[0] = \
                    int(current_children[0].name) - 1
            else:
                node_list[num_nodes - (distances[node.name] + 1)].children[0] = \
                    num_nodes - (distances[current_children[0].name] + 1)

            # Child 1
            if not current_children[1].children:
                # if current_children[1].is_leaf():
                node_list[num_nodes - (distances[node.name] + 1)].children[1] = \
                    int(current_children[1].name) - 1
            else:
                node_list[num_nodes - (distances[node.name] + 1)].children[1] = \
                    num_nodes - (distances[current_children[1].name] + 1)
    return TREE(num_leaves, node_list, -1)


def findpath_distance(t1, t2, c=False):
    # This being called with two ete3 trees
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]

    if c:
        return lib.findpath_distance(t1, t2)

    # ct1 = tree_io.read_newick(t1.write(format=5))
    # ct2 = tree_io.read_newick(t2.write(format=5))
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)

    return lib.findpath_distance(ct1, ct2)


# TODO decorator for findpath_distance to check type of given arguments and return the
#  correct function


class TREE_LIST(Structure):
    _fields_ = [('num_trees', c_int), ('trees', POINTER(TREE))]

    # def __init_(self, num_trees, trees):
    #     self.num_trees = num_trees
    #     self.trees = trees


# C function tree_to_string (for testing purposes)
tree_to_cluster_string = lib.tree_to_string
tree_to_cluster_string.argtypes = [POINTER(TREE)]
tree_to_cluster_string.restype = c_char_p

# C function return_findpath
findpath_path = lib.return_findpath
findpath_path.argtypes = [POINTER(TREE), POINTER(TREE)]
findpath_path.restype = TREE_LIST
