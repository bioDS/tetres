# algorithm for condiditonal clade distrbitution by Bret Larget
from collections import defaultdict
from decimal import *


def get_maps(trees):

    m1 = defaultdict(int)  # map for each clade how often it got sampled
    m2 = defaultdict(int)  # map for each (c1,c2) clade how often this specific relation got sampled
    uniques = {}

    seen = {}
    # todo m1 needs to have the root and the single leaf clades
    #  need to get a list of weights for the unique trees
    #  this needs to replace trees in the for loop and the weight will be added instead of 1 each time!

    for ix, t in enumerate(trees):
        # if t.get_newick(f=9) == '(8,(((((1,10),7),2),3),((5,4),(6,9))));' or t.get_newick(f=9) == '(8,(((((1,10),7),2),3),((5,4),(9,6))));':
        #     print("Galah")
        if not frozenset(sorted(t.get_clades())) in seen:
            seen[frozenset(sorted(t.get_clades()))] = ix
            uniques[ix] = [] 
            # uniques.append(ix)
        else:
            # todo get the index of where current tree has been seen before!
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
