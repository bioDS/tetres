from tetres.partitions.conditional_partitions import *


def test_get_conditional_partitions(ten_taxa_multichain):
    get_conditional_partitions(ten_taxa_multichain[0].trees[0])
    assert True


def test_get_list_of_partitions(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    assert True


def test_get_pp(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    pp = get_pp(ten_taxa_multichain[0].trees[0], dict_part)
    assert True


def test_get_pp_coverage_threespace(threespace):
    cov = get_pp_coverage(threespace[0].trees)
    assert cov == 1, "Coverage failed (3space)"


def test_get_pp_coverage_fourspace(fourspace):
    cov = get_pp_coverage(fourspace.trees)
    cov == 1, "Coverage failed (4space)!"


def test_get_greedy_pp_tree(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    tree = get_greedy_pp_tree(dict_part, len(ten_taxa_multichain[0].trees[0]))
    assert True


def test_are_compatible(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    keys = list(dict_part.keys())
    assert are_compatible(keys[1], keys[2]), "Are compatible function failed"


def test_are_compatible2(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    keys = list(dict_part.keys())
    assert not are_compatible(keys[1], keys[3]), "Are compatible function failed"


def test_get_tree_from_partition(ten_taxa_multichain):
    distance_sum = 0
    for tree in ten_taxa_multichain[0].trees:
        cp = get_conditional_partitions(tree)
        cp = [cp[k] for k in sorted(cp, reverse=True)[1:]]  # adapting for now
        nt = get_tree_from_partition(cp, 10)
        distance_sum += nt.fp_distance(tree)
    assert distance_sum == 0, "Tree from partition failed!"


def test_get_tree_from_partition20(twenty_taxa_tts):
    cp = get_conditional_partitions(twenty_taxa_tts[0])
    cp = [cp[k] for k in sorted(cp, reverse=True)[1:]]  # adapting for now
    nt = get_tree_from_partition(cp, len(twenty_taxa_tts[0]))
    assert nt.fp_distance(twenty_taxa_tts[0]) == 0, "Get tree form partition for 20 taxa failed!"


def test_sample_from_dict_partition(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    t = sample_from_dict_partition(dict_part, samples=100, n_taxa=len(ten_taxa_multichain[0].trees[0]))
    assert True


def test_search_maxpp_tree(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    t, t_lp = search_maxpp_tree(dict_part, n_taxa=len(ten_taxa_multichain[0].trees[0]))
    baseline = get_greedy_pp_tree(dict_part, len(ten_taxa_multichain[0].trees[0]))
    baseline_p = get_pp(baseline, dict_part, log=True)
    assert t_lp > baseline_p, "Branch and bound search went wrong."


def test_get_greedy_relaxed_pp_tree(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    relaxed_tree = get_greedy_relaxed_pp_tree(dict_part, len(ten_taxa_multichain[0].trees[0]))
    greedy_tree = get_greedy_pp_tree(dict_part, len(ten_taxa_multichain[0].trees[0]))
    assert relaxed_tree.fp_distance(greedy_tree) == 1, "Something went wrong with relaxed greedy!"


def test_convert_to_relaxed_partition_dict(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    relaxed_dict = convert_to_relaxed_partition_dict(dict_part)
    greedy_tree = get_greedy_pp_tree(dict_part, len(ten_taxa_multichain[0].trees[0]))
    relaxed_greedy_tree = get_greedy_pp_tree(relaxed_dict, len(ten_taxa_multichain[0].trees[0]))
    relaxed_tree = get_greedy_relaxed_pp_tree(dict_part, len(ten_taxa_multichain[0].trees[0]))
    assert (greedy_tree.fp_distance(relaxed_tree), greedy_tree.fp_distance(relaxed_greedy_tree), relaxed_greedy_tree.fp_distance(relaxed_tree)) == (1, 1, 0), "Converting dict of partitions to relaxed dict of partitions failed!"


def test_calc_Entropy_cpd(twenty_taxa_tts):
    dict_part = get_dict_of_partitions(twenty_taxa_tts)
    h = calc_Entropy_cpd(dict_part)
    assert h ==  16.116557383124142, "Entropy Calculation failed!"


def test_max_cpd_tree_dp(twenty_taxa_tts):
    dict_part = get_dict_of_partitions(twenty_taxa_tts)
    # greedy = get_greedy_pp_tree(dict_part, 20)
    dp_tree, dp_prob = max_cpd_tree_dp(dict_part)
    test_p = get_pp(dp_tree, dict_part)
    assert test_p == dp_prob, "DP alg for partitions failed!"


def test_max_cpd_tree_dp_greedy(twenty_taxa_tts):
    dict_part = get_dict_of_partitions(twenty_taxa_tts)
    greedy = get_greedy_pp_tree(dict_part, 20)
    dp_tree, dp_prob = max_cpd_tree_dp(dict_part)
    greedy_p = get_pp(dp_tree, dict_part)
    assert greedy_p <= dp_prob, "DP alg for partitions failed!"


def test_add_tree_dict_partitions(ten_taxa_multichain):
    dict_part = get_dict_of_partitions(ten_taxa_multichain[0].trees)
    add_tree(dict_part, ten_taxa_multichain[0].trees[0])
    assert dict_part[sorted(dict_part.keys(), key=len)[0]]["Count"] == 1002, "Failed to add tree to dictionary of partitions"


def test_relaxed_dp(twenty_taxa_tts):
    dict_part = get_dict_of_partitions(twenty_taxa_tts)
    count, relaxed_dict = convert_to_relaxed_partition_dict(dict_part, ret_newpaths_count=True)
    dp_tree, dp_prob = max_cpd_tree_dp(relaxed_dict)
    assert dp_prob > 0, "Failed Relaxed DP."
