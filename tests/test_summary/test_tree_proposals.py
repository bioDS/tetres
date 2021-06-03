
from treeoclock.summary._tree_proposals import search_neighbourhood_greedy
from treeoclock.trees.time_trees import TimeTreeSet


def test_search_neighbourhood_greedy_ncorenone(five_taxa_tts):
    res = search_neighbourhood_greedy(t=five_taxa_tts[1], trees=five_taxa_tts, t_value=10)
    assert (five_taxa_tts[0].fp_distance(res[0]), res[1]) == (0, 5), "Greedy neighbourhood search with n_cores=None failed!"


def test_search_neighbourhood_greedy_ncore1(five_taxa_tts):
    res = search_neighbourhood_greedy(t=five_taxa_tts[1], trees=five_taxa_tts, t_value=10, n_cores=1)
    assert (five_taxa_tts[0].fp_distance(res[0]), res[1]) == (0, 5), "Greedy neighbourhood search with n_cores=1 failed!"


def test_search_neighbourhood_greedy_ncore4(five_taxa_tts):
    res = search_neighbourhood_greedy(t=five_taxa_tts[1], trees=five_taxa_tts, t_value=10, n_cores=4)
    assert (five_taxa_tts[0].fp_distance(res[0]), res[1]) == (0, 5), "Greedy neighbourhood search with n_cores=4 failed!"

# TODO test for the selection methods, first construct a test case with multiple choices and then test those
#  also test for the value error that gets raised if wrong select is given
