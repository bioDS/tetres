from tetres.online.online_distance import greedy_online_omp, search_neighbourhood_greedy_online_omp, intelligent_neighbourhood_online

import tracemalloc
import gc


def test_intelligent_neighbourhood_online_twelve(twelve_taxa_tts):
    twelve_taxa_tts.get_common_clades()
    # todo add memory usage check for the test
    # starting the monitoring
    tracemalloc.start()
    mem_before = tracemalloc.get_tracemalloc_memory()
    for t in twelve_taxa_tts[1:100]:
        cur_move_list = intelligent_neighbourhood_online(t, twelve_taxa_tts.common_clades)
        assert cur_move_list.num_moves <= 2 * (len(t)-2), "Too many moves in move_list!"
    gc.collect()
    mem_after = tracemalloc.get_tracemalloc_memory()
    tracemalloc.stop()
    print(mem_after, mem_before)
    assert mem_after - mem_before < 50000, "Failed Memory usage test (12 taxa)"


def test_intelligent_neighbourhood_online_twenty(twenty_taxa_tts):
    twenty_taxa_tts.get_common_clades()
    # todo add memory usage check for the test
    # starting the monitoring
    tracemalloc.start()
    mem_before = tracemalloc.get_tracemalloc_memory()
    for t in twenty_taxa_tts[1:100]:
        cur_move_list = intelligent_neighbourhood_online(t, twenty_taxa_tts.common_clades)
        assert cur_move_list.num_moves <= 2 * (len(t) - 2), "Too many moves in move_list!"
    gc.collect()
    mem_after = tracemalloc.get_tracemalloc_memory()
    tracemalloc.stop()
    assert mem_after - mem_before < 50000, "Failed Memory usage test (20 taxa)"




def test_greedy_online_omp(twelve_taxa_tts, twelve_taxa_tts_start):
    temp_ret = greedy_online_omp(trees=twelve_taxa_tts, start=twelve_taxa_tts_start[0])
    assert True, "Failed!"


# todo i need a bunch of tests for the new funcitons now
#  check whether it finds the same best neighbour
#  check what happens if we have no better neighbour and how to deal with that
#  how do we store the centroid, just moves or do we iteratively update it? not sure
#  then I will also need to refactor the c code where these test will be very helpful

