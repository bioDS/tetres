from treeoclock.trees._converter import ctree_to_ete3
from treeoclock.trees.time_trees import TimeTree
import re


def get_conditional_partitions(tree, as_dict=False):
    rank_etree = ctree_to_ete3(tree.ctree)
    rc_dict = {}
    for node in rank_etree.traverse("levelorder"):
        if len(node) > 1:
            c = node.children
            c0_leafs = []
            for leaf in c[0]:
                c0_leafs.append(int(leaf.name))
            c1_leafs = []
            for leaf in c[1]:
                c1_leafs.append(int(leaf.name))
            node_rank = int(node.name.replace("I", ""))

            # c0 should always be the one with lowest integer among the two lists
            if min(c1_leafs) < min(c0_leafs):
                c1_leafs, c0_leafs = c0_leafs, c1_leafs

            rc_dict[node_rank] = f"{','.join([str(i) for i in sorted(c0_leafs)])}|{','.join([str(i) for i in sorted(c1_leafs)])}"

    rc_dict[max(rc_dict.keys())+1] = ",".join(sorted([t for t in re.split(r",|\|", rc_dict[max(rc_dict.keys())])], key=int))  # Setting root

    if as_dict:
        return rc_dict

    for rank in sorted(rc_dict.keys(), reverse=True)[1:]:
        cur_split = rc_dict[rank]
        to_be_mod = rc_dict[rank+1]
        to_be_swapped = ','.join(sorted([t for t in re.split(r',|\|', cur_split)], key=int))
        rc_dict[rank] = re.sub(to_be_swapped, cur_split, to_be_mod)
    # return set of strings for a tree
    return rc_dict


def get_dict_of_partitions(treeset):
    ret = {}
    for t in treeset:
        t_part = get_conditional_partitions(t)
        for rank in range(max(t_part.keys()), 2, -1):
            # t_part[rank+1]
            if t_part[rank] in ret:
                if t_part[rank - 1] in ret[t_part[rank]]:
                    ret[t_part[rank]][t_part[rank - 1]] += 1
                    ret[t_part[rank]]["Count"] += 1
                else:
                    ret[t_part[rank]][t_part[rank - 1]] = 1
                    ret[t_part[rank]]["Count"] += 1
            else:
                ret[t_part[rank]] = {t_part[rank - 1]: 1, "Count": 1}
    return ret


def get_pp(tree, dict_partitions):
    tree_part = get_conditional_partitions(tree)
    pp = 1
    for k in sorted(tree_part.keys(), reverse=True)[:-2]:
        if tree_part[k] in dict_partitions:
            if tree_part[k-1] in dict_partitions[tree_part[k]]:
                pp *= (dict_partitions[tree_part[k]][tree_part[k-1]]/dict_partitions[tree_part[k]]["Count"])
            else:
                return 0
        else:
            return 0
    return pp


def get_pp_coverage(trees):
    dict_partitions = get_dict_of_partitions(trees)
    cov = 0
    # todo this should only count unique trees and nothing else
    for t in trees:
        cov += get_pp(t, dict_partitions)
    return cov


# todo general sample from dict_partition function


def get_greedy_pp_tree(dict_partitions, n_taxa):
    out = []
    # todo this is dependant on the fact that the dict is sorted, should be sorted by number of | in string
    k = sorted(dict_partitions.keys(), key=len)[0]
    while len(out) < n_taxa-2:
        highest = 0
        h_v = 0
        for v in dict_partitions[k]:
            if v != "Count":
                if dict_partitions[k][v] > highest:
                    highest = dict_partitions[k][v]
                    h_v = v
        k = h_v
        out.append(h_v)
    # calculate a tree from the list of partitions
    t = get_tree_from_partition(out, n_taxa)
    return out


import ete3


def get_tree_from_partition(p_list, n_taxa):
    cur_t = ete3.Tree(support=n_taxa-1, name=",".join([str(i) for i in range(1, n_taxa+1)]))
    # add the first thing manually
    init_split = p_list[0].split("|")
    for split in init_split:
        if len(split) == 1:
            # It splits off a leaf, therefore the distance can be set to the rank n-1
            cur_t.add_child(name=split, dist=n_taxa-1)
        else:
            cur_t.add_child(name=split)
    rank = n_taxa-2
    for string in p_list[1:]:
        splits = string.split("|")
        splits = [s for s in splits if s not in cur_t]

        # find the node we need to extend:
        node = cur_t.search_nodes(name=",".join(sorted(",".join(splits).split(","), key=int)))[0]
        node.support = rank
        node.dist = node.up.support - rank
        for s in splits:
            if len(s) == 1:
                node.add_child(name=s, dist=rank)
            else:
                node.add_child(name=s)
        rank -= 1
    return TimeTree(cur_t.write(format=0))
