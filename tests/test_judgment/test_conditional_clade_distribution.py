from treeoclock.judgment.conditional_clade_distribution import get_maps, get_tree_probability
from collections import Counter
from decimal import *


def test_get_maps(ten_taxa_cMChain):
    m1, m2 = get_maps(ten_taxa_cMChain[0].trees)
    assert True # todo not sure


def test_get_tree_probability(ten_taxa_cMChain):
    m1, m2 = get_maps(ten_taxa_cMChain[0].trees)
    p = get_tree_probability(ten_taxa_cMChain[0].trees[1], m1, m2)
    assert p == 0.0037693679383000523, "Get tree probability failed!"


def test_coverage(ten_taxa_cMChain):
    m1, m2, uniques = get_maps(ten_taxa_cMChain[0].trees)
    probs = []
    for u in uniques.keys():
        probs.append(get_tree_probability(ten_taxa_cMChain[0].trees[u], m1, m2))
    # The sum over the probabilities of all unique trees should be 1 (if all of the distribution or the whole space has been sampled)
    assert sum(Counter(probs)) == 0.8403013465493053, "Coverage sum failed"


def test_coverage2(twenty_taxa_tts):
    m1, m2 = get_maps(twenty_taxa_tts)
    probs = []
    for t in twenty_taxa_tts:
        probs.append(get_tree_probability(t, m1, m2))
    # The sum over the probabilities of all unique trees should be 1 (if all of the distribution or the whole space has been sampled)
    assert sum(Counter(probs)) == 0.6749181338892793, "Coverage sum failed"

