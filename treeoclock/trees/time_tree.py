import re
import os
import ete3
from ctypes import POINTER, CDLL, c_long

from treeoclock.trees._converter import ete3_to_ctree, ctree_to_ete3
from treeoclock.trees._ctrees import TREE
from treeoclock.trees.findpath_distance import findpath_distance

# TODO temporary imports
# import line_profiler

# TODO doumentation
#  How/Where to put away all the functions so that this file only contains the two classes ?
#  Maybe even split up the two classes in two files ?!
#  And then all the timetree handling options are done via the classes and its functions
from treeoclock.trees.time_tree_set import TimeTreeSet

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')


class TimeTree:
    def __init__(self, nwk):
        self.etree = ete3.Tree(nwk)
        self.ctree = ete3_to_ctree(self.etree)

    def __len__(self):
        return self.ctree.num_leaves

    def fp_distance(self, tree):
        return findpath_distance(self.ctree, tree.ctree)

    def get_newick(self, f=5):
        return self.etree.write(format=f)

    def write_newick(self, file_path, f=5):
        # TODO maybe make this a write_nexus function with the correct mapping dict and all ?
        return self.etree.write(outfile=file_path, format=f)

    # TODO one_neighbourhood function


# TODO
#  Write trees with TimeTreeSet.write(index=) or index range or just write() to get all trees as nexus output ?
#  len() function for TimeTreeSet
#  list of trees should be numpy.array(dtype=TimeTree)
#  Function for remapping, i.e. changing the mapping dict and changing each tree in the list


def neighbourhood(tree):
    lib.rank_move.argtypes = [POINTER(TREE), c_long]

    num_leaves = len(tree)

    # All possible rank moves
    print(tree.get_newick(f=9))
    for i in range(num_leaves - 1):
        lib.rank_move(tree.ctree, i)  # TODO Changes the tree inline, but that should not be the case!
        print(ctree_to_ete3(tree.ctree).write(format=9))

    # All NNI moves
    lib.nni_move.argtypes = [POINTER(TREE), c_long]

    # TODO how to check if nni move is possible for a given index i and i+1 ?
    #  Maybe get exceptions from the c code ? or do the same logic here so that the c code does not fail
    # This should return a numpy array with TimeTree objects
    return 0


if __name__ == '__main__':
    d_name = 'Dengue'

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    # TODO testing of findpath_path and also the TimeTree.fp_path() function ?!
    # TODO do some timing comparison with direct calling, i.e. no type recognition, much difference ?
    print(myts[1].fp_distance(myts[0]))
    print(findpath_distance(myts[1], myts[0]))
    print(findpath_distance(myts[1].ctree, myts[0].ctree))
    print(findpath_distance(myts[1].etree, myts[0].etree))
    
    # n = neighbourhood(myts[0])

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
