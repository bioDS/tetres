__author__ = 'Lena Collienne and Jordan Kettles'

import os
from ctypes import *  # TODO

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


def ete3_to_ctree(tree):
    num_leaves = len(tree.get_tree_root())
    num_nodes = (num_leaves * 2) - 1
    node_list = (NODE * num_nodes)()

    distances = []

    tree_root = tree.get_tree_root()

    for node in tree.traverse('levelorder'):
        if not node.is_leaf():
            distances.append(node.get_distance(tree_root))

    distances = sorted(distances)
    if not len(set(distances)) == num_leaves-1:
        sys.exit('Distances to root not unique! \n'
                 'This has to be resolved!')
    for node in tree.traverse('levelorder'):
        if node.is_leaf():
            node_list[int(node.name) - 1].parent = num_nodes - (distances.index(node.up.get_distance(tree_root)) + 1)
        else:
            if node.is_root():
                node_list[num_nodes - (distances.index(node.get_distance(tree_root)) + 1)].parent = num_nodes - 1
            else:
                node_list[num_nodes - (distances.index(node.get_distance(tree_root)) + 1)].parent = \
                    num_nodes - (distances.index(node.up.get_distance(tree_root)) + 1)
            current_children = node.get_children()
            if len(current_children) != 2:
                sys.exit('Not a binary tree!')

            # Child 0
            if current_children[0].is_leaf():
                node_list[num_nodes - (distances.index(node.get_distance(tree_root)) + 1)].children[0] = \
                    int(current_children[0].name) - 1
            else:
                node_list[num_nodes - (distances.index(node.get_distance(tree_root)) + 1)].children[0] = \
                    num_nodes - (distances.index(current_children[0].get_distance(tree_root)) + 1)

            # Child 1
            if current_children[1].is_leaf():
                node_list[num_nodes - (distances.index(node.get_distance(tree_root)) + 1)].children[1] = \
                    int(current_children[1].name) - 1
            else:
                node_list[num_nodes - (distances.index(node.get_distance(tree_root)) + 1)].children[1] = \
                    num_nodes - (distances.index(current_children[1].get_distance(tree_root)) + 1)
    return TREE(num_leaves, node_list, -1)


def findpath_distance(t1, t2, c=False):
    # This being called with two ete3 trees
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]

    import tree_io

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
