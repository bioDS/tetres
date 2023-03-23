import pytest

from tetres.summary._tree_proposals import search_neighbourhood_greedy, NoBetterNeighbourFound, search_neighbourhood_greedy_omp
from tetres.trees.time_trees import TimeTree


def test_search_neighbourhood_greedy_ncorenone(five_taxa_tts):
    res = search_neighbourhood_greedy(t=five_taxa_tts[1], trees=five_taxa_tts, t_value=10)
    assert (five_taxa_tts[0].fp_distance(res[0]), res[1]) == (0, 5), "Greedy neighbourhood search with n_cores=None failed!"


def test_search_neighbourhood_greedy_ncore1(five_taxa_tts):
    res = search_neighbourhood_greedy(t=five_taxa_tts[1], trees=five_taxa_tts, t_value=10, n_cores=1)
    assert (five_taxa_tts[0].fp_distance(res[0]), res[1]) == (0, 5), "Greedy neighbourhood search with n_cores=1 failed!"


def test_search_neighbourhood_greedy_ncore4(five_taxa_tts):
    res = search_neighbourhood_greedy(t=five_taxa_tts[1], trees=five_taxa_tts, t_value=10, n_cores=4)
    assert (five_taxa_tts[0].fp_distance(res[0]), res[1]) == (0, 5), "Greedy neighbourhood search with n_cores=4 failed!"


def test_search_neighbourhood_greedy_nothingfound_exception(five_taxa_tts):
    with pytest.raises(NoBetterNeighbourFound):
        # This Error is raised because of t_value being to low!
        other = TimeTree('((1:3,5:3):1,(3:2,(4:1,2:1):1):2);')
        search_neighbourhood_greedy(t=other, trees=five_taxa_tts, t_value=10)


def test_search_neighbourhood_greedy_valueerror_selection(five_taxa_tts):
    with pytest.raises(ValueError):
        search_neighbourhood_greedy(t=five_taxa_tts[0], trees=five_taxa_tts, t_value=10, select='Test')


def test_search_neighbourhood_greedy_openmp(five_taxa_tts):
    res = search_neighbourhood_greedy_omp(t=five_taxa_tts[1], trees=five_taxa_tts, t_value=10, n_cores=1)
    assert (five_taxa_tts[0].fp_distance(res[0]), res[1]) == (0, 5), "Greedy neighbourhood search with openmp failed!"
