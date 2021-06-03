from random import shuffle
from typing import Union

from treeoclock.trees.time_trees import TimeTreeSet, findpath_path, TimeTree
from treeoclock.trees._converter import ctree_to_ete3


def frechet_mean(trees: Union[TimeTreeSet, list]):

    working_tree_set = trees.copy()  # Making a copy of the list of trees
    shuffle(working_tree_set)  # Shuffle the copy
    frechet_mean = working_tree_set.pop().ctree  # Initialize the FM tree with the first tree in the list

    fm_divider = 2  # Starting with the tree at the middle of the path i.e. 1/2 of the length

    while working_tree_set:
        # While there are still trees left in the set
        cur_tree = working_tree_set.pop()
        # cur_path = frechet_mean.fp_path(cur_tree)
        cur_path = findpath_path(frechet_mean, cur_tree.ctree)
        # Only consider to break if the tree isn't close to the current FM tree
        if int(len(cur_path) / fm_divider) == 0:
            break
        frechet_mean = cur_path[int(len(cur_path)/fm_divider)]
        # Assigning the new FM as the tree 1/fm_divider of the distances to the cur_tree
        fm_divider += 1

    return TimeTree(ctree_to_ete3(frechet_mean).write(format=5))


if __name__ == '__main__':
    d_name = "RSV2"

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    for _ in range(5):
        fm = frechet_mean(myts)
        # print(fm.get_newick())

    for _ in range(5):
        fm = frechet_mean(myts[1:2])
        # print(fm.get_newick())
