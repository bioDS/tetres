from tetres.trees._converter import ctree_to_ete3
from tetres.trees.time_trees import TimeTree
import re
import ete3
import numpy as np


def get_unconditional_partitions(tree):

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

    # currently no need for the full set of leaves in this dict
    # rc_dict[max(rc_dict.keys())+1] = ",".join(sorted([t for t in re.split(r",|\|", rc_dict[max(rc_dict.keys())])], key=int))  # Setting root

    # todo fix bug that the values do not contain all the taxa!!!

    for rank in sorted(rc_dict.keys(), reverse=True)[2:]:
        cur_split = rc_dict[rank]
        to_be_mod = rc_dict[rank+1]
        to_be_swapped = ','.join(sorted([t for t in re.split(r',|\|', cur_split)], key=int))
        rc_dict[rank] = re.sub(fr"{to_be_swapped}(?=$|[^\d])",cur_split, to_be_mod)
    # return set of strings for a tree
    return rc_dict


def get_dict_of_unconditional_partitions(treeset):
    ret = {}
    for t in treeset:
        t_part = get_unconditional_partitions(t)
        for rank in t_part:
            if rank in ret:
                if t_part[rank] in ret[rank]:
                    ret[rank][t_part[rank]] += 1
                else:
                    ret[rank][t_part[rank]] = 1
            else:
                ret[rank] = {t_part[rank]: 1}
    return ret


def get_uncond_pp(tree, dict_partitions, tree_size, log=False):
    tree_part = get_unconditional_partitions(tree)
    if log:
        pp=0
    else:
        pp = 1
    for k in tree_part:
        if tree_part[k] in dict_partitions[k]:
            if log:
                pp += np.log(dict_partitions[k][tree_part[k]]/tree_size)
            else:
                pp *= (dict_partitions[k][tree_part[k]]/tree_size)
        else:
            return None
    return pp


def get_uncond_pp_coverage(trees):
    dict_partitions = get_dict_of_unconditional_partitions(trees)
    cov = 0
    # todo this should only count unique trees and nothing else
    for t in trees:
        cov += get_uncond_pp(t, dict_partitions, tree_size=len(trees))
    return cov


def get_greedy_uncond_pp_tree(dict_partitions, n_taxa):
    out = []
    for k in dict_partitions:
        cur_max = 0
        mv = 0
        for v in dict_partitions[k]:
            cur_val = dict_partitions[k][v]
            if cur_val > cur_max:
                cur_max = cur_val
                mv = v
        out.append(mv)
    # append last h_v with replaced , for | to get full resolution building the tree
    out.append(re.sub(",", "|", h_v))
    # calculate a tree from the list of partitions
    return(get_tree_from_partition(out, n_taxa))


def get_tree_from_partition(p_list, n_taxa):

    # todo probably need to adapt something here, like add the lowest and highest rank partitions
    #  but those are trivial anyway ....

    cur_t = ete3.Tree(support=n_taxa-1, name=",".join([str(i) for i in range(1, n_taxa+1)]))
    # add the first thing manually
    init_split = p_list[0].split("|")
    rank = n_taxa-1
    for split in init_split:
        if len(split.split(",")) == 1:
            # It splits off a leaf, therefore the distance can be set to the rank n-1
            cur_t.add_child(name=split, dist=rank)
        else:
            cur_t.add_child(name=split)
    rank -= 1
    for string in p_list[1:]:
        splits = string.split("|")
        splits = [s for s in splits if s not in cur_t]
        # find the node we need to extend:
        node = cur_t.search_nodes(name=",".join(sorted(",".join(splits).split(","), key=int)))[0]
        node.support = rank
        node.dist = node.up.support - rank
        for s in splits:
            if len(s.split(",")) == 1:
                node.add_child(name=s, dist=rank)
            else:
                node.add_child(name=s)
        rank -= 1
    return TimeTree(cur_t.write(format=0))


def sample_from_dict_partition(dict_partitions, n_taxa, samples=1):
    # sample samples many trees from a given dict_partition
    new_samples = []
    for _ in range(samples):
        out = []
        # todo this is dependent on the fact that the dict is sorted, should be sorted by number of | in string
        k = sorted(dict_partitions.keys(), key=len)[0]
        while len(out) < n_taxa - 2:
            cur_dict = dict_partitions[k].copy()
            cur_count = cur_dict.pop("Count")  # delete count from dict and get the total number of choices
            cur_p = [i / cur_count for i in list(cur_dict.values())]
            out.append(np.random.choice(list(cur_dict.keys()), p=cur_p))
            k = out[-1]
        out.append(re.sub(",", "|", out[-1]))  # todo this probably needs to be sorted?
        new_samples.append(get_tree_from_partition(out, n_taxa))
    return new_samples


def search_maxpp_tree(dict_partitions, n_taxa):
    baseline = get_greedy_pp_tree(dict_partitions, n_taxa)
    baseline_p = get_pp(baseline, dict_partitions, log=True)

    sol = {}

    def rec_max_search(log_p, partition, level):
        nonlocal baseline_p
        nonlocal dict_partitions
        nonlocal sol
        if log_p < baseline_p:
            return None
        if level == 2:
            sol[log_p] = []
            return log_p
        for k in dict_partitions[partition]:
            if k != "Count":
                lp = rec_max_search(log_p+np.log(dict_partitions[partition][k]/dict_partitions[partition]["Count"]), k, level-1)
                if lp != None:
                    # print(dict_partitions[partition][k], dict_partitions[partition]["Count"])
                    sol[lp].append(k)
                    return lp
                else:
                    return None
    out = rec_max_search(0, list(dict_partitions.keys())[0], level=n_taxa)
    if out == None:
        return None, None
    # todo check if there is actually anything in sol otherwise return -1 or something
    sol[max(sol)].insert(0, sol[max(sol)][0].replace(",", "|"))
    sol[max(sol)].reverse()
    t = get_tree_from_partition(sol[max(sol)], n_taxa)
    return t, max(sol)