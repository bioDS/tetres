from treeoclock.trees._ctrees import TREE_LIST, TREE
from treeoclock.trees.time_trees import TimeTreeSet, findpath_path, TimeTree
from treeoclock.trees._converter import ctree_to_ete3
from treeoclock.summary.compute_sos import compute_sos_mt
from multiprocessing.pool import ThreadPool as Pool

from random import shuffle
from typing import Union
import os
from ctypes import CDLL, POINTER
import warnings

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
        with warnings.catch_warnings():
            # Ignores the 'Free memory' warning issued by findpath_path
            warnings.simplefilter("ignore")
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


def fm_sos(trees: TimeTreeSet):
    fm = frechet_mean(trees)
    sos = compute_sos_mt(fm, trees)
    return fm, sos


def loop_fm(trees: TimeTreeSet, n_fms=10, n_cores:int = None):
    with Pool(n_cores) as p:
        fms = p.map(fm_sos, [trees for _ in range(n_fms)])
    return min(fms, key=lambda t: t[1])[0]


# # todo think about a timeTreeset without c trees or a list of ete3 trees and then convert to ctrees when needed
# def multiple_frechet_mean(tree_file: str, n_fms=4, n_cores=4):
#     with mp.Pool(n_cores) as p:
#         fms = p.map(fm_sos, [tree_file for _ in range(n_fms)])
#     return TimeTree(min(fms, key=lambda t: t[1])[0])


if __name__ == '__main__':
    d_name = 'RSV2'

    tree_file = f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees'
    myts = TimeTreeSet(tree_file)

    fm = loop_fm(myts)

    print(fm.get_newick())
