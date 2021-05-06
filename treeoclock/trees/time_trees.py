import re
import os
import ete3
import functools
from ctypes import POINTER, CDLL, c_long

from _converter import ete3_to_ctree, ctree_to_ete3
from _ctrees import TREE, TREE_LIST

# TODO temporary imports
# import line_profiler

# TODO doumentation
#  How/Where to put away all the functions so that this file only contains the two classes ?
#  Maybe even split up the two classes in two files ?!
#  And then all the timetree handling options are done via the classes and its functions

lib = CDLL(f'{os.path.dirname(os.path.realpath(__file__))}/findpath.so')


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


class TimeTreeSet:
    def __init__(self, file):
        self.map = get_mapping_dict(file)
        self.trees = my_trees_read(file)

    def __getitem__(self, index):
        return self.trees[index]

    def __len__(self):
        return len(self.trees)


# TODO
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


@functools.singledispatch
def findpath_distance(arg, *args):
    raise TypeError(type(arg) + " not supported.")


@findpath_distance.register(TREE)
def findpath_distance_c(t1, t2):
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    return lib.findpath_distance(t1, t2)


@findpath_distance.register(ete3.Tree)
def findpath_distance_ete3(t1, t2):
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    ct1 = ete3_to_ctree(t1)
    ct2 = ete3_to_ctree(t2)
    return lib.findpath_distance(ct1, ct2)


@findpath_distance.register(TimeTree)
def findpath_distance_ete3(t1, t2):
    lib.findpath_distance.argtypes = [POINTER(TREE), POINTER(TREE)]
    return lib.findpath_distance(t1.ctree, t2.ctree)


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
    d_name = 'Dengue'

    # t = read_nexus(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees', c=False)
    # ct = read_nexus(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees', c=True)

    # myt = my_trees_read(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    print(myts[1].fp_distance(myts[0]))
    print(findpath_distance(myts[1], myts[0]))
    print(findpath_distance(myts[1].ctree, myts[0].ctree))
    print(findpath_distance(myts[1].etree, myts[0].etree))

    findpath_distance(1, 2)
    
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
