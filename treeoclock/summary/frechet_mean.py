from treeoclock.trees._ctrees import TREE_LIST, TREE
from treeoclock.trees.time_trees import TimeTreeSet, findpath_path, TimeTree
from treeoclock.trees._converter import ctree_to_ete3

from random import shuffle
from typing import Union
import os
from ctypes import CDLL, POINTER

lib = CDLL(f"{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/trees/findpath.so")


def frechet_mean(trees: Union[TimeTreeSet, list]):
    lib.copy_tree.argtypes = [POINTER(TREE)]
    lib.copy_tree.restype = TREE  # Deep copy tree
    lib.free_tree.argtypes = [POINTER(TREE)]  # free tree memory
    lib.free_treelist.argtypes = [TREE_LIST]  # free tree_list memory

    index_list = list(range(0, len(trees)))  # indices of all trees

    shuffle(index_list)  # Shuffle the copy
    frechet_mean = lib.copy_tree(trees[index_list.pop()].ctree)  # Initialize the FM tree with the first tree in the list

    fm_divider = 2  # Starting with the tree at the middle of the path i.e. 1/2 of the length

    while index_list:
        # While there are still trees left in the set
        cur_tree = trees[index_list.pop()]
        path = findpath_path(frechet_mean, cur_tree.ctree)

        # Only consider to break if the tree isn't close to the current FM tree
        pos = int(path.num_trees / fm_divider)
        if pos == 0:
            break
        lib.free_tree(frechet_mean)  # Free the memory of the previous Frechet Mean tree as it is a deep copy
        # Assigning the new FM as the tree 1/fm_divider of the distances to the cur_tree
        frechet_mean = lib.copy_tree(path.trees[pos])
        lib.free_treelist(path)  # Free the returned path Tree_List
        fm_divider += 1
    ret_tree = TimeTree(ctree_to_ete3(frechet_mean).write(format=5))  # Make a TimeTree to return
    lib.free_tree(frechet_mean)  # Free the memory of the Frechet mean tree
    return ret_tree


if __name__ == '__main__':
    d_name = 'RSV2'

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    fm, mem = frechet_mean(myts)
    # print(fm.get_newick())

    import matplotlib.pyplot as plt
    plt.plot(mem)
    plt.show()
