# algorithm for condiditonal clade distrbitution by Bret Larget
from collections import defaultdict
from decimal import *
import ete3
import math
from tetres.trees.time_trees import TimeTree


def get_maps(trees):

    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled
    uniques = {}

    seen = {}
    # todo m1 needs to have the root and the single leaf clades
    #  need to get a list of weights for the unique trees
    #  this needs to replace trees in the for loop and the weight will be added instead of 1 each time!

    for ix, t in enumerate(trees):
        if not frozenset(sorted(t.get_clades())) in seen:
            seen[frozenset(sorted(t.get_clades()))] = ix
            uniques[ix] = []
        else:
            uniques[seen[frozenset(sorted(t.get_clades()))]].append(ix)

        for node in t.etree.traverse("levelorder"):
            if len(node) > 2:
                c = node.children
                c0_leafs = set()
                for leaf in c[0]:
                    c0_leafs.add(int(leaf.name))
                c1_leafs = set()
                for leaf in c[1]:
                    c1_leafs.add(int(leaf.name))
                parent_clade = frozenset(sorted(c0_leafs.union(c1_leafs)))
                m1[parent_clade] += 1
                if min(c0_leafs) < min(c1_leafs):
                    m2[(parent_clade, frozenset(c0_leafs))] += 1
                else:
                    m2[(parent_clade, frozenset(c1_leafs))] += 1
    return m1, m2, uniques


def get_tree_probability(tree, m1, m2):
    # getcontext().prec = 20
    probability = 1
    for node in tree.etree.traverse("levelorder"):
        if len(node) > 2:
            c = node.children
            c0_leafs = set()
            for leaf in c[0]:
                c0_leafs.add(int(leaf.name))
            c1_leafs = set()
            for leaf in c[1]:
                c1_leafs.add(int(leaf.name))
            parent_clade = frozenset(sorted(c0_leafs.union(c1_leafs)))
            if m1[parent_clade] != 0:
                if min(c0_leafs) < min(c1_leafs):
                    probability *= Decimal(m2[(parent_clade, frozenset(c0_leafs))]) / Decimal(m1[parent_clade])
                else:
                    probability *= Decimal(m2[(parent_clade, frozenset(c1_leafs))]) / Decimal(m1[parent_clade])
        if probability == 0:
            return 0
    return float(probability)


def get_ccd_tree_branch_bound(m1, m2, prob):

    # todo figure out why this returns the same tree multiple times
    #  todo make some improvements to the efficiency, maybe add multiprocessing???

    working_list = [([max(m1.keys())], [], 1)]  # initialize with root clade, empty tree, probability 1
    output = []
    while working_list:
        # Get the currently unresolved tree as cur_tuple
        cur_tuple = working_list.pop()
        for tuple_split in cur_tuple[0]:
            # iterate over all unresolved parts of the tree

            all_splits_of_cur_tuple = [(m2[p, c], c) for p, c in [(list(i)[0], list(i)[1]) for i in m2.keys()] if
                                        p == tuple_split]
            # todo here we only need to look at splits once, because if we look at it a second time then it will end with the same result?!?!?!

            for split in all_splits_of_cur_tuple:
                # iterate over all possible splits of this currently unresolved part of the tree
                new_prob = cur_tuple[2] * (split[0] / m1[tuple_split])
                if new_prob >= prob:
                    # if the probability of this resolving is still higher than the one given as threhsold
                    # add it to the working list of partially unresolved trees
                    new_list = cur_tuple[0].copy()  # copy all unresolved bits of tree
                    new_list.remove(tuple_split)  # remove bit that is resolved in this step
                    if len(split[1]) > 2:
                        new_list.append(split[1])  # add child1 if it contains more than 2 taxe, i.e. is unresolved
                    if len(tuple_split.difference(split[1])) > 2:
                        new_list.append(tuple_split.difference(split[1]))  # add child2 if it is unresolved
                    tree_as_list_of_splits = cur_tuple[1].copy()  # Copy the tree which is less resolved
                    tree_as_list_of_splits.append((tuple_split, split[1]))  # add this pair parent, child1 split to the output tree list
                    if not new_list:
                        # if the new list does not contain anything we have a fully resolved tree
                        #  hence we add it to output and don't need to resolve it further

                        # todo check that output does not contain the current element, i.e. make output a list?
                        output.append((tree_as_list_of_splits,
                                         new_prob))

                    else:
                        # current tree is not yet fully resolved, so add it to working list

                        # todo same thing as above, check whether this new list is not already in there
                        #  i.e. is the same form of tree already being worked on, if yes, abandon this branch instead
                        # if tree_as_list_of_splits not in [i[1] for i in working_list]:
                        working_list.append((new_list,
                                             tree_as_list_of_splits,
                                             new_prob))
    best_bb_tree = max(output, key=lambda item: item[1])
    return get_tree_from_list_of_splits(best_bb_tree[0]), output, best_bb_tree[1]


def get_ccd_tree_dfs(m1, m2, prob):
    # working_list = [([max(m1.keys())], [], 1)]  # initialize with root clade, empty tree, probability 1

    useless = {}  # saving all clades that will not lead to a better tree
    seen_resolved_clades = {}

    all_clades = sorted(list(m1.keys()), key=len)
    for current_clade in all_clades:
        all_splits_current = [i for i in m2 if i[0] == current_clade]
        for split in all_splits_current:
            # this for loop needs to find the best split of the current parent, if any exists!
            child1 = split[1]
            child2 = split[0].difference(split[1])

            if child1 in useless or child2 in useless:
                # Current Split is not going to improve solution
                continue

            if len(child1) < 3:
                # length 2 or 1 gives probability 1
                c1_prob = 1
            else:
                if child1 in seen_resolved_clades:
                    c1_prob, _ = seen_resolved_clades[child1]
            if len(child2) < 3:
                # length 2 or 1 gives probability 1
                c2_prob = 1
            else:
                if child2 in seen_resolved_clades:
                    c2_prob, _ = seen_resolved_clades[child2]

            if (c1_prob * c2_prob) < prob:
                # The current split will not lead to an improvement
                continue

            cur_prob = m2[split] / m1[split[0]]  # Prob of current parent, given split
            split_prob = c1_prob * c2_prob * cur_prob  # best probability of current parent with split
            if math.isclose(split_prob, prob) or split_prob > prob:
            # if split_prob >= prob:  # does not work due to float comparison, 0.1 + 0.2 != 0.3, use math.isclose instead of ==
                if split[0] in seen_resolved_clades:
                    # parent was already found, do we need to update?
                    if seen_resolved_clades[split[0]][0] <= split_prob:
                        # this split has better probability, therefore update the seen_resolved_clades with the better split of split[0]
                        seen_resolved_clades[split[0]] = (split_prob, child1)
                else:
                    # we have not seen the parent before
                    seen_resolved_clades[split[0]] = (split_prob, child1)
        if current_clade not in seen_resolved_clades:
            # all possible splits of current_clade have not resulted in a tree with higher probability
            #  therefore add it to the useless pile
            useless[current_clade] = 0
    # todo build a tree from the seen and resolved dict
    return seen_resolved_clades


def get_greedy_ccd_tree(m1, m2):
    working_list = [max(m1.keys())]  # initialize with root clade
    greedy_tree = []
    while working_list:
        cur_parent = working_list.pop()
        all_splits_of_cur_parent = [(m2[p,c], c) for p, c in [(list(i)[0], list(i)[1]) for i in m2.keys()] if p == cur_parent]
        # todo maybe add instant probability calculation at some point too
        cur_greedy_split = max(all_splits_of_cur_parent, key=lambda  item: item[0])[1]
        if len(cur_greedy_split) > 2:
            working_list.append(cur_greedy_split)
        greedy_tree.append((cur_parent, cur_greedy_split))
        if len(cur_parent.difference(cur_greedy_split)) > 2:
            # if the second child is not a resolved tree yet (i.e. 2 or 1 taxa) add to working list
            working_list.append(cur_parent.difference(cur_greedy_split))
    return greedy_tree


def get_tree_from_list_of_splits(splits):
    # this is dependent on the current structure of how the greedy list of splits is created!
    n_taxa = len(splits[0][0])
    dist = 1
    # support is the rank of the node...
    cur_t = ete3.Tree(support=0, dist=0, name=",".join([str(i) for i in sorted(splits[0][0])]))
    for parent, child1 in splits:
        node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(parent)]))[0]
        child2 = parent.difference(child1)
        if len(child1) == 1:
            node.add_child(name=",".join([str(i) for i in sorted(child1)]), support=n_taxa-1)
        else:
            node.add_child(name= ",".join([str(i) for i in sorted(child1)]), support=dist)
            dist +=1
        if len(child2) == 1:
            node.add_child(name=",".join([str(i) for i in sorted(child2)]), support=n_taxa-1)
        else:
            node.add_child(name= ",".join([str(i) for i in sorted(child2)]), support=dist)
            dist +=1
        # IF the children nodes have 2 taxa the correspoinding leafs need to be added here
        if len(child1) == 2:
            node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(child1)]))[0]
            node.add_child(name=list(sorted(child1))[0], support=n_taxa-1)
            node.add_child(name=list(sorted(child1))[1], support=n_taxa-1)
        if len(child2) == 2:
            node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(child2)]))[0]
            node.add_child(name=list(sorted(child2))[0], support=n_taxa-1)
            node.add_child(name=list(sorted(child2))[1], support=n_taxa-1)

    # setting the node distances to ranks, no meaning just to have a ranked tree
    for node in cur_t.iter_descendants("postorder"):
        node.dist = node.support - node.up.support

    return TimeTree(cur_t.write(format=5))


def add_centroid(centroid, m1, m2):
    for node in centroid.etree.traverse("levelorder"):
        if len(node) > 2:
            c = node.children
            c0_leafs = set()
            for leaf in c[0]:
                c0_leafs.add(int(leaf.name))
            c1_leafs = set()
            for leaf in c[1]:
                c1_leafs.add(int(leaf.name))
            parent_clade = frozenset(sorted(c0_leafs.union(c1_leafs)))
            if m1[parent_clade] == 0:
                m1[parent_clade] += 10
                if min(c0_leafs) < min(c1_leafs):
                    m2[(parent_clade, frozenset(c0_leafs))] += 1
                else:
                    m2[(parent_clade, frozenset(c1_leafs))] += 1
            else:
                x = int(m1[parent_clade]/9)
                m1[parent_clade] += x
                if min(c0_leafs) < min(c1_leafs):
                    m2[(parent_clade, frozenset(c0_leafs))] += x
                else:
                    m2[(parent_clade, frozenset(c1_leafs))] += x

    return m1, m2
