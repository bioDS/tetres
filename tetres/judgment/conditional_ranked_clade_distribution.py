# algorithm for condiditonal clade distrbitution by Bret Larget
from collections import defaultdict
from decimal import *
from tetres.trees._converter import ctree_to_ete3


def get_rank_maps(trees):

    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled
    uniques = {}

    seen = {}
    # todo m1 needs to have the root and the single leaf clades
    #  need to get a list of weights for the unique trees
    #  this needs to replace trees in the for loop and the weight will be added instead of 1 each time!

    for ix, t in enumerate(trees):
        rank_etree = ctree_to_ete3(t.ctree)
        # todo all unique rank trees?
        # if not frozenset(sorted(t.get_clades())) in seen:
        #     seen[frozenset(sorted(t.get_clades()))] = ix
        #     uniques[ix] = []
        # else:
        #     uniques[seen[frozenset(sorted(t.get_clades()))]].append(ix)

        for node in rank_etree.traverse("levelorder"):
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
                if len(c0_leafs) == 1:
                    if len(c1_leafs) == 1:
                        # both leafs, we do the regular
                        if min(c0_leafs) < min(c1_leafs):
                            m2[(parent_clade, frozenset(c0_leafs))] += 1
                        else:
                            m2[(parent_clade, frozenset(c1_leafs))] += 1
                    else:
                        # one leaf, the other vertex has higher rank
                        m2[(parent_clade, frozenset(c1_leafs))] += 1
                elif len(c1_leafs) == 1:
                    # c1 is leaf but c0 is not!
                    m2[(parent_clade, frozenset(c0_leafs))] += 1
                else:
                    # both no leafs, extract name and get one with higher rank
                    c0_rank = int(c[0].name.replace("I", ""))
                    c1_rank = int(c[1].name.replace("I", ""))
                    if c0_rank > c1_rank:
                        m2[(parent_clade, frozenset(c0_leafs))] += 1
                    else:
                        m2[(parent_clade, frozenset(c1_leafs))] += 1

                # if min(c0_leafs) < min(c1_leafs):
                #     m2[(parent_clade, frozenset(c0_leafs))] += 1
                # else:
                #     m2[(parent_clade, frozenset(c1_leafs))] += 1

    return m1, m2, uniques


def get_rank_tree_probability(tree, m1, m2):
    # getcontext().prec = 20
    probability = 1
    rank_etree = ctree_to_ete3(tree.ctree)
    for node in rank_etree.traverse("levelorder"):
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
                if len(c0_leafs) == 1:
                    if len(c1_leafs) == 1:
                        # both leafs, we do the regular
                        if min(c0_leafs) < min(c1_leafs):
                            probability *= Decimal(m2[(parent_clade, frozenset(c0_leafs))]) / Decimal(m1[parent_clade])
                        else:
                            probability *= Decimal(m2[(parent_clade, frozenset(c1_leafs))]) / Decimal(m1[parent_clade])
                    else:
                        # one leaf, the other vertex has higher rank
                        probability *= Decimal(m2[(parent_clade, frozenset(c1_leafs))]) / Decimal(m1[parent_clade])
                elif len(c1_leafs) == 1:
                    # c1 is leaf but c0 is not!
                    probability *= Decimal(m2[(parent_clade, frozenset(c0_leafs))]) / Decimal(m1[parent_clade])
                else:
                    # both no leafs, extract name and get one with higher rank
                    c0_rank = int(c[0].name.replace("I", ""))
                    c1_rank = int(c[1].name.replace("I", ""))
                    if c0_rank > c1_rank:
                        probability *= Decimal(m2[(parent_clade, frozenset(c0_leafs))]) / Decimal(m1[parent_clade])
                    else:
                        probability *= Decimal(m2[(parent_clade, frozenset(c1_leafs))]) / Decimal(m1[parent_clade])
        if probability == 0:
            return 0
    return float(probability)


def add_rank_centroid(centroid, m1, m2):
    raise ValueError("Not sure about this implementation anymore! Why is it added 10 times? does that explain the high probability?")
    rank_centroid_etree = ctree_to_ete3(centroid.ctree)
    for node in rank_centroid_etree.traverse("levelorder"):
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
