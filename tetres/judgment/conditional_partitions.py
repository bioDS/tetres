from tetres.trees._converter import ctree_to_ete3
from tetres.trees.time_trees import TimeTree
import re
import ete3
import numpy as np
import copy


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
        rc_dict[rank] = re.sub(fr"{to_be_swapped}(?=$|[^\d])",cur_split, to_be_mod)
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


def get_pp(tree, dict_partitions, log=False):
    tree_part = get_conditional_partitions(tree)
    if log:
        pp=0
    else:
        pp = 1
    for k in sorted(tree_part.keys(), reverse=True)[:-2]:
        if tree_part[k] in dict_partitions:
            if tree_part[k-1] in dict_partitions[tree_part[k]]:
                if log:
                    pp += np.log(dict_partitions[tree_part[k]][tree_part[k-1]]/dict_partitions[tree_part[k]]["Count"])
                else:
                    pp *= (dict_partitions[tree_part[k]][tree_part[k-1]]/dict_partitions[tree_part[k]]["Count"])
            else:
                return 0
        else:
            return 0
    return pp


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
    # append last h_v with replaced , for | to get full resolution building the tree
    out.append(re.sub(",", "|", h_v))
    # calculate a tree from the list of partitions
    return(get_tree_from_partition(out, n_taxa))


def are_compatible(partition1, partition2):
    set1 = set(partition1.split("|"))
    set2 = set(partition2.split("|"))
    d1 = set1 -set2
    d2 = set2 - set1
    if len(d1) > len(d2):
        # switching so that d1 is the smaller set
        d1, d2 = d2, d1
    if len(d1) !=1:
        return False
    if len(d2) != 2:
        return False
    d2_union = d2.pop() + "," + d2.pop()
    return ",".join(sorted(d2_union.split(','), key=int)) in d1


def convert_to_relaxed_partition_dict(dict_partitions, ret_newpaths_count=False):
    relaxed_dict = copy.deepcopy(dict_partitions)
    new_paths = 0
    # todo implement this and test it vs the greedy relaxed one ....
    for fix_partition in dict_partitions:
        for other_partition in dict_partitions:
            if fix_partition != other_partition:
                if len(fix_partition.split("|")) == len(other_partition.split("|")):
                    for sub_partition in dict_partitions[other_partition]:
                        if sub_partition != "Count":
                            if are_compatible(sub_partition, fix_partition):
                                new_paths += 1
                                if sub_partition in relaxed_dict[fix_partition]:
                                    relaxed_dict[fix_partition][sub_partition] += dict_partitions[other_partition][sub_partition]
                                else:
                                    relaxed_dict[fix_partition][sub_partition] = dict_partitions[other_partition][sub_partition]
                                relaxed_dict[fix_partition]["Count"] += dict_partitions[other_partition][sub_partition]
    # print(f"NewPaths: {new_paths}")
    if ret_newpaths_count:
        return new_paths, relaxed_dict
    return relaxed_dict


def get_greedy_relaxed_pp_tree(dict_partitions, n_taxa):
    out = []
    # todo this is dependant on the fact that the dict is sorted, should be sorted by number of | in string
    out.append(sorted(dict_partitions.keys(), key=len)[0])
    while len(out) < n_taxa-1:
        cur_dict = dict_partitions[out[-1]].copy()
        cur_dict.pop("Count")
        # find all rank 2 keys in dict_partitions
        for other_k in dict_partitions:
            if other_k != "Count" and other_k != out[-1]:
                if len(other_k.split("|")) == len(out[-1].split("|")):
                    for sub_other_k in dict_partitions[other_k]:
                        if are_compatible(sub_other_k, out[-1]):
                            if sub_other_k in cur_dict:
                                cur_dict[sub_other_k] += dict_partitions[other_k][sub_other_k]
                            else:
                                cur_dict[sub_other_k] = dict_partitions[other_k][sub_other_k]
        # todo this is not relevant here, but I need this for the calculation of probabilities so the datastructure needs to be adapted for these relaxed probabilities
        #  probably just calculate the dicitonary like this in the first place so that i can use the same greedy choice algorithm than before?
        highest = 0
        h_v = 0
        for k in cur_dict:
            if cur_dict[k] > highest:
                    highest = cur_dict[k]
                    h_v = k
        out.append(h_v)
    # append last h_v with replaced , for | to get full resolution building the tree
    out.append(re.sub(",", "|", out[-1]))
    # calculate a tree from the list of partitions
    return get_tree_from_partition(out[1:], n_taxa)


def get_tree_from_partition(p_list, n_taxa):
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
    if out is None:
        return None, None
    sol[max(sol)].insert(0, sol[max(sol)][0].replace(",", "|"))
    sol[max(sol)].reverse()
    t = get_tree_from_partition(sol[max(sol)], n_taxa)
    return t, max(sol)