import re

from treeoclock.trees.time_tree import TimeTree
from treeoclock.trees.findpath_distance import findpath_distance
from treeoclock.trees.findpath_path import findpath_path


class TimeTreeSet:
    def __init__(self, file):
        self.map = get_mapping_dict(file)
        self.trees = read_nexus(file)

    def __getitem__(self, index):
        return self.trees[index]

    def __len__(self):
        return len(self.trees)
    
    def findpath_distance(self, i, j):
        return findpath_distance(self.trees[i].ctree, self.trees[j].ctree)
    
    def findpath_path(self, i, j):
        return findpath_path(self.trees[i].ctree, self.trees[j].ctree)


# TODO
#  Write trees with TimeTreeSet.write(index=) or index range or just write() to get all trees as nexus output ?
#  Function for remapping, i.e. changing the mapping dict and changing each tree in the list


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