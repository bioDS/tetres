from tetres.online.online_distance import greedy_online_omp, search_neighbourhood_greedy_online_omp,\
    intelligent_neighbourhood_online, get_precomp_list, free_precomp_list, _compute_and_free_precomp_list
from tetres.summary._tree_proposals import search_neighbourhood_greedy, NoBetterNeighbourFound
from tetres.summary.compute_sos import compute_sos_mt
from tetres.trees.time_trees import TimeTree

import pytest
import psutil
import gc


def test_intelligent_neighbourhood_online_twelve(twelve_taxa_tts):
    twelve_taxa_tts.get_common_clades()
    before = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    for t in twelve_taxa_tts[1:100]:
        cur_move_list = intelligent_neighbourhood_online(t, twelve_taxa_tts.common_clades)
        assert cur_move_list.num_moves <= 2 * (len(t)-2), "Too many moves in move_list!"
    after = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    assert after - before < 10.0, "Failed Memory usage test (12 taxa)"


def test_intelligent_neighbourhood_online_twenty(twenty_taxa_tts):
    twenty_taxa_tts.get_common_clades()
    before = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    for t in twenty_taxa_tts[1:100]:
        cur_move_list = intelligent_neighbourhood_online(t, twenty_taxa_tts.common_clades)
        assert cur_move_list.num_moves <= 2 * (len(t) - 2), "Too many moves in move_list!"
    after = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    assert after - before < 10.0, "Failed Memory usage test (20 taxa)"


def test_suspicious_memory(twelve_taxa_tts, twelve_taxa_tts_start):
    import tracemalloc
    before = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    tracemalloc.start()
    s0 = tracemalloc.take_snapshot()
    for _ in range(10):
        _compute_and_free_precomp_list(twelve_taxa_tts, twelve_taxa_tts_start[0])
    after = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    s1 = tracemalloc.take_snapshot()
    assert after-before < 10.0, "Failed memory usage test"


def test_get_precomp_list_twelve_memory(twelve_taxa_tts, twelve_taxa_tts_start):
    before = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    for _ in range(100):
        precomp_list = get_precomp_list(twelve_taxa_tts, twelve_taxa_tts_start[0])
        free_precomp_list(precomp_list, len(twelve_taxa_tts_start[0]))
    gc.collect()
    after = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    # 1995.28515625, with freeing precomps, still memory leak, not sure where
    # 2798.62890625, without freeing precomps
    assert after-before < 10.0, "get_precomp_list failed Memory usage Test (12 taxa)"


def test_get_precomp_list_twenty(twenty_taxa_tts, twenty_taxa_tts_start):
    before = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    for _ in range(100):
        precomp_list = get_precomp_list(twenty_taxa_tts, twenty_taxa_tts_start[0])
    after = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    assert after - before < 10.0, "get_precomp_list failed Memory usage Test (12 taxa)"


def test_memory_search_neighbourhood_greedy_online_omp_twelve(twelve_taxa_tts, twelve_taxa_tts_start):
    precomp_list = get_precomp_list(twelve_taxa_tts, twelve_taxa_tts_start[0])
    before = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    centroid, sos = search_neighbourhood_greedy_online_omp(twelve_taxa_tts_start[0], twelve_taxa_tts, precomp_list)
    after = psutil.Process().memory_info().rss / 1024 ** 2  # in MiB
    assert after-before < 10.0, "search_neighbourhood_greedy_online_omp failed Memory usage! There is a memory leak."


def test_functionality_search_neighbourhood_greedy_online_omp_twelve(twelve_taxa_tts, twelve_taxa_tts_start):
    precomp_list = get_precomp_list(twelve_taxa_tts, twelve_taxa_tts_start[0])
    best_neigh_online, sos_online = search_neighbourhood_greedy_online_omp(twelve_taxa_tts_start[0], twelve_taxa_tts, precomp_list)
    start_sos = compute_sos_mt(twelve_taxa_tts_start[0], twelve_taxa_tts)
    best_neigh_offline, sos_offline = search_neighbourhood_greedy(twelve_taxa_tts_start[0], twelve_taxa_tts, t_value=start_sos)
    assert all([sos_online == sos_offline, best_neigh_online.fp_distance(best_neigh_offline) == 0]),\
        "Failed to find the same best neighbour using the online distance for greedy neighbourhood search!"


# todo test what if we do not have a better neighbour, implement throwing the exception
def test_nobetterneighbour_search_neighbourhood_greedy_online_omp_twelve(five_taxa_tts):
    with pytest.raises(NoBetterNeighbourFound):
        # This Error is raised because of t_value being to low!
        other = TimeTree('((1:3,5:3):1,(3:2,(4:1,2:1):1):2);')
        precomp_list = get_precomp_list(five_taxa_tts, other)
        search_neighbourhood_greedy_online_omp(other, five_taxa_tts, precomp_list)


def test_search_neighbourhood_greedy_online_omp_twenty(twelve_taxa_tts, twenty_taxa_tts_start):

    assert True


def test_greedy_online_omp(twelve_taxa_tts, twelve_taxa_tts_start):
    temp_ret = greedy_online_omp(trees=twelve_taxa_tts, start=twelve_taxa_tts_start[0])
    assert True, "Failed!"


# todo i need a bunch of tests for the new funcitons now
#  check whether it finds the same best neighbour
#  check what happens if we have no better neighbour and how to deal with that
#  how do we store the centroid, just moves or do we iteratively update it? not sure
#  then I will also need to refactor the c code where these test will be very helpful

