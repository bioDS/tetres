import re
import os
import ete3
from ctypes import POINTER, CDLL

from _converter import ete3_to_ctree
from _ctrees import TREE, TREE_LIST
# TODO temporary imports
import line_profiler

# TODO doumentation for the classes missing

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')


class TimeTree:
    def __init__(self, nwk):
        self.etree = ete3.Tree(nwk)
        self.ctree = ete3_to_ctree(self.etree)

    def fp_distance(self, tree):
        return findpath_distance(self.ctree, tree.ctree, True)

    def get_newick(self):
        return self.etree.write(format=5)

    def write_newick(self, file):
        # TODO maybe make this a write_nexus function with the correct mapping dict and all ?
        return self.etree.write(outfile=file, format=5)

    def neighbourhood(self):

        # This should return a numpy array with TimeTree objects
        return 0

    # TODO one_neighbourhood function


class TimeTreeSet:
    def __init__(self, file):
        self.map = get_mapping_dict(file)
        self.trees = read_nexus(file)

# TODO class should be iterable
#  the initialization with a nexus, reading the trees
#  and saving the mapping dictionary within
#  access the list of trees with TimeTreeSet.trees ?
#  access the mapping dict with TimeTreeSet.map
#  Write trees with TimeTreeSet.write(index=) or index range or just write() to get all trees as nexus output ?
#  len() function for TimeTreeSet
#  list of trees should be numpy.array(dtype=TimeTree)
#  Function for remapping, i.e. changing the mapping dict and changing each tree in the list


def findpath_path(t1, t2, c=False):
    # TODO maybe a decorator to get rid of the c parameter ?
    # C function return_findpath, returns python list of TREE objects
    lib.return_findpath.argtypes = [POINTER(TREE), POINTER(TREE)]
    lib.return_findpath.restype = TREE_LIST
    if c:
        path = lib.return_findpath(t1, t2)
        return [path.trees[i] for i in range(path.num_trees)]
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    path = lib.return_findpath(ct1, ct2)
    return [path.trees[i] for i in range(path.num_trees)]


def findpath_distance(t1, t2, c=False):
    # TODO maybe a decorator to get rid of the c parameter ?
    # This being called with two ete3 trees
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]

    if c:
        return lib.findpath_distance(t1, t2)
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    return lib.findpath_distance(ct1, ct2)


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

    import sys
    import random
    from timeit import default_timer as timer

    d_name = 'Dengue'

    t = read_nexus(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees', c=False)
    
    ct = read_nexus(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees', c=True)
    
    myt = my_trees_read(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')
    
    from _converter import ctree_to_ete3
    
    ctree_to_ete3(ct[0])
    
    p = findpath_path(ct[0], ct[1], True)
    ep = findpath_path(t[0], t[1], False)

    # from pympler import asizeof
    #
    # print(asizeof.asizeof(t))
    # print(asizeof.asizeof(ct))
    # print(asizeof.asizeof(myt))
