# algorithm for condiditonal clade distrbitution by Bret Larget
from collections import defaultdict
from decimal import *
import ete3
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


def get_greedy_ccd_tree(m1, m2):
    # todo somehow return a greedy tree from the maps of partitions

    working_list = [max(m1.keys())]
    greedy_tree = []
    while working_list:
        cur_parent = working_list.pop()
        all_splits_of_cur_parent = [(m2[p,c], c) for p, c in [(list(i)[0], list(i)[1]) for i in m2.keys()] if p == cur_parent]
        # todo maybe add instant probability calculation at some point too
        cur_greedy_split = max(all_splits_of_cur_parent, key=lambda  item: item[0])[1]
        if len(cur_greedy_split) > 2:
            working_list.append(cur_greedy_split)
            # todo add the other child, i.e. cur_parent-cur_greedy_split
        greedy_tree.append((cur_parent, cur_greedy_split))
        if len(cur_parent.difference(cur_greedy_split)) > 2:
            # if the second child is not a resolved tree yet (i.e. 2 or 1 taxa) add to working list
            working_list.append(cur_parent.difference(cur_greedy_split))
    return greedy_tree


def get_tree_from_list_of_splits(splits):
    print("HELP ME")
    # this is dependent on the current structure of how the greedy list of splits is created!
    n_taxa = len(splits[0])
    # support is the rank of the node...
    cur_t = ete3.Tree(support=n_taxa - 1, name=",".join([str(i) for i in sorted(splits[0][0])]))
    for parent, child1 in splits:
        node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(parent)]))[0]
        child2 = parent.difference(child1)
        node.add_child(name= ",".join([str(i) for i in sorted(child1)]))
        node.add_child(name= ",".join([str(i) for i in sorted(child2)]))
        # IF the children nodes have 2 taxa the correspoinding leafs need to be added here
        if len(child1) == 2:
            node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(child1)]))[0]
            node.add_child(name=list(sorted(child1))[0])
            node.add_child(name=list(sorted(child1))[1])
        if len(child2) == 2:
            node = cur_t.search_nodes(name=",".join([str(i) for i in sorted(child2)]))[0]
            node.add_child(name=list(sorted(child2))[0])
            node.add_child(name=list(sorted(child2))[1])
    # todo annotate with branch lengths here?
    # todo otherwise how do we get it to be a timetree? or do we even need that at all?
    # return TimeTree(cur_t.write(format=0))
    return cur_t.write(format=0)


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
