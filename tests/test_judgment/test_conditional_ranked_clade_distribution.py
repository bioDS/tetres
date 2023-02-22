from treeoclock.judgment.conditional_ranked_clade_distribution import get_rank_maps, get_rank_tree_probability
from collections import Counter


def test_get_maps(ten_taxa_cMChain):
    m1, m2, uniques = get_rank_maps(ten_taxa_cMChain[0].trees)
    assert True # todo not sure


def test_get_tree_probability(ten_taxa_cMChain):
    m1, m2, uniques = get_rank_maps(ten_taxa_cMChain[0].trees)
    p = get_rank_tree_probability(ten_taxa_cMChain[0].trees[1], m1, m2)
    assert p == 0.0029781866133667644, "Get tree probability failed!"


def test_coverage(ten_taxa_cMChain):
    m1, m2, uniques = get_rank_maps(ten_taxa_cMChain[0].trees)
    probs = []
    # for u in uniques.keys():
    for u in ten_taxa_cMChain[0].trees:
        probs.append(get_rank_tree_probability(u, m1, m2))
    # The sum over the probabilities of all unique trees should be 1 (if all of the distribution or the whole space has been sampled)
    # todo only sum over unique trees
    assert sum(probs) == 0.6759011979706053, "Coverage sum failed"


def test_coverage2(twenty_taxa_tts):
    m1, m2, uniques = get_rank_maps(twenty_taxa_tts)
    probs = []
    for t in twenty_taxa_tts:
        probs.append(get_rank_tree_probability(t, m1, m2))
    # The sum over the probabilities of all unique trees should be 1 (if all of the distribution or the whole space has been sampled)
    # todo only sum over unique trees
    assert sum(probs) == 0.3937814210255327, "Coverage sum failed"

def test_coverage3(fourspace):
    m1, m2, uniques = get_rank_maps(fourspace[0].trees)
    probs = []
    for t in fourspace[0].trees:
        probs.append(get_rank_tree_probability(t, m1, m2))
    # The sum over the probabilities of all unique trees should be 1 (if all of the distribution or the whole space has been sampled)
    assert sum(probs) == 1, "Coverage sum failed"


def test_coverage3space(threespace):
    m1, m2, uniques = get_rank_maps(threespace[0].trees)
    probs = []
    for t in threespace[0].trees:
        probs.append(get_rank_tree_probability(t, m1, m2))
    # The sum over the probabilities of all unique trees should be 1 (if all of the distribution or the whole space has been sampled)
    assert sum(probs) == 1, "Coverage sum failed"
