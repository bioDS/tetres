from tetres.online.online_distance import greedy_online_omp, search_neighbourhood_greedy_online_omp,\
    intelligent_neighbourhood_online, get_precomp_list, free_precomp_list

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


def test_search_neighbourhood_greedy_online_omp_twelve(twelve_taxa_tts, twelve_taxa_tts_start):
    search_neighbourhood_greedy_online_omp(twelve_taxa_tts_start, twelve_taxa_tts)

    assert True


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

