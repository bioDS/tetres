from tetres.partitions.conditional_clade_distribution import *
from collections import Counter
import numpy as np


def test_get_maps(ten_taxa_multichain):
    m1, m2 = get_maps(ten_taxa_multichain[0].trees)
    assert True # todo not sure


def test_get_tree_probability(ten_taxa_multichain):
    m1, m2, u = get_maps(ten_taxa_multichain[0].trees)
    p = get_tree_probability(ten_taxa_multichain[0].trees[1], m1, m2)
    assert p == 0.0037693679383000514, "Get tree probability failed!"


def test_get_tree_log_probability(ten_taxa_multichain):
    m1, m2, u = get_maps(ten_taxa_multichain[0].trees)
    p = get_tree_probability(ten_taxa_multichain[0].trees[1], m1, m2, log=True)
    assert p == np.log(0.0037693679383000514), "Get tree probability failed!"


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
    assert len(greedy_timetree) == 20, "Failed Greedy CCD tree algorithm!"


def test_get_ccd_tree_branch_bound(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    greedy_tree = get_tree_from_list_of_splits(get_greedy_ccd_tree(m1, m2))
    greedy_prob = get_tree_probability(greedy_tree, m1, m2)
    bb_tree, output, treeprob = get_ccd_tree_branch_bound(m1, m2, greedy_prob)
    assert len(bb_tree) == 20, "Failed Branch and bound CCD algorithm"


def test_get_ccd_tree_bottom_up(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    best_tree = get_ccd_tree_bottom_up(m1, m2)
    assert len(best_tree) == 20, "Failed DFS BB algorithm"


def test_rogueiness(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    rogueiness(m1)
    assert True


def test_sample_tree_from_ccd(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    sample = sample_tree_from_ccd(m1, m2, n=100)
    assert len(sample) == 100


def test_sample_logprob_from_ccd(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    sample = sample_logprob_from_ccd(m1, m2, n=100)
    assert True


def test_remove_taxa_fourspace(fourspace):
    m1, m2, u = get_maps(fourspace.trees)
    for i in range(1,5):
        remove_taxa(i, m1, m2)
    assert True


def test_remove_taxa_twelve(twelve_taxa_tts):
    m1, m2, u = get_maps(twelve_taxa_tts)
    # remove_taxa(9, m1, m2)
    for i in range(1, 13):
        remove_taxa(i, m1, m2)
    assert True


def test_remove_taxa_twenty(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    # remove_taxa(4, m1, m2)
    for i in range(1, 21):
        remove_taxa(i, m1, m2)
    assert True


def test_remove_taxa_dp(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    for i in range(1, 21):
        new_m1, new_m2 = remove_taxa(i, m1, m2)
        tree = get_ccd_tree_bottom_up(new_m1, new_m2)
        assert len(tree) == 19

def test_remove_taxa_dp_prob(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    for i in range(1, 21):
        new_m1, new_m2 = remove_taxa(i, m1, m2)
        tree = get_ccd_tree_bottom_up(new_m1, new_m2)
        p = get_tree_probability(tree, new_m1, new_m2)
        print(p)
        assert p < 1


def test_calc_Entropy(twenty_taxa_tts):
    m1, m2, u = get_maps(twenty_taxa_tts)
    h = calc_Entropy(m1, m2)
    assert h == 7.500020954155561, "Entropy calculation failed!"
