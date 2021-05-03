__author__ = 'Lena Collienne, Lars Berling, Jordan Kettles'

import re
import sys
import ete3
from collections import OrderedDict
from call_findpath import TREE, NODE, TREE_LIST, findpath_distance
from ctypes import c_long

# TODO temporary imports
import line_profiler


def get_mapping_dict(file: str) -> dict:
    """
    Returns the taxon mapping of the nexus file as a dictionary

    :param file: A nexus file path
    :type file: str
    :return: Dictionary containing the mapping of taxa(values) to int(keys)
    :rtype: dict {int --> str}
    """
    # Extract a mapping dict from a file dict: int --> taxa
    # Begin trees;
    #   Translate

    # ;
    # trees
    # End;
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
                if not c:
                    tree_string = f'{re.split(re_tree, line)[1][:re.split(re_tree, line)[1].rfind(")")+1]};'
                    trees.append(ete3.Tree(re.sub(brackets, "", tree_string)))
                else:
                    tree_string = f'{re.split(re_tree, line)[1][:re.split(re_tree, line)[1].rfind(")")+1]};'
                    trees.append(read_newick(re.sub(brackets, "", tree_string)))
    return trees


if __name__ == '__main__':

    # TODO need one_neighbourhood
    # TODO need conversion from cTREE to a ete3 tree
    # TODO TREELIST to ete3 trees
    # TODO TREELIST extract one tree with an index ?

    import sys
    from timeit import default_timer as timer
    d_name = 'RSV2'

    from call_findpath import ctree_to_ete3, ete3_to_ctree

    t = read_nexus(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees', c=False)
    cct = [ete3_to_ctree(tr) for tr in t]
    etec = [ctree_to_ete3(tr) for tr in cct]




