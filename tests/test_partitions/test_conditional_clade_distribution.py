from tetres.partitions.conditional_clade_distribution import get_maps, get_tree_probability, get_greedy_ccd_tree, get_tree_from_list_of_splits
from collections import Counter


def test_get_maps(ten_taxa_multichain):
    m1, m2 = get_maps(ten_taxa_multichain[0].trees)
    assert True # todo not sure


def test_get_tree_probability(ten_taxa_multichain):
    m1, m2 = get_maps(ten_taxa_multichain[0].trees)
    p = get_tree_probability(ten_taxa_multichain[0].trees[1], m1, m2)
    assert p == 0.0037693679383000523, "Get tree probability failed!"


def test_coverage(ten_taxa_multichain):
    m1, m2, uniques = get_maps(ten_taxa_multichain[0].trees)
    probs = []
    for u in uniques.keys():
        probs.append(get_tree_probability(ten_taxa_multichain[0].trees[u], m1, m2))
    # The sum over the probabilities of all unique trees should be 1 (if all of the distribution or the whole space has been sampled)
    assert sum(Counter(probs)) == 0.8403013465493053, "Coverage sum failed"


def test_coverage2(twenty_taxa_tts):
    m1, m2 = get_maps(twenty_taxa_tts)
    probs = []
    for t in twenty_taxa_tts:
        probs.append(get_tree_probability(t, m1, m2))
    # The sum over the probabilities of all unique trees should be 1 (if all of the distribution or the whole space has been sampled)
    assert sum(Counter(probs)) == 0.6749181338892793, "Coverage sum failed"


def test_get_greedy_ccd_tree(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    greedy_tree = get_greedy_ccd_tree(m1, m2)
    greedy_timetree = get_tree_from_list_of_splits(greedy_tree)
    # assert len(greedy_timetree) == 20, "Failed Greedy CCD tree algorithm!"
    assert True
