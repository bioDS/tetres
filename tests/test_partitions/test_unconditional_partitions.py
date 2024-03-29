from tetres.partitions.unconditional_partitions import get_unconditional_partitions, get_dict_of_unconditional_partitions, get_uncond_pp,get_uncond_pp_coverage, get_greedy_uncond_pp_tree


def test_get_unconditional_partitions(ten_taxa_multichain):
    rc_dict = get_unconditional_partitions(ten_taxa_multichain[0].trees[0])
    assert True


def test_get_list_of_unconditional_partitions(ten_taxa_multichain):
    dict_part = get_dict_of_unconditional_partitions(ten_taxa_multichain[0].trees)
    assert True


def test_get_uncond_pp(ten_taxa_multichain):
    dict_part = get_dict_of_unconditional_partitions(ten_taxa_multichain[0].trees)
    pp = get_uncond_pp(ten_taxa_multichain[0].trees[100], dict_part, tree_size=len(ten_taxa_multichain[0].trees), log=False)
    log_pp = get_uncond_pp(ten_taxa_multichain[0].trees[100], dict_part, tree_size=len(ten_taxa_multichain[0].trees), log=True)
    assert True


# todo this is a weird case for the current implementation
# def test_get_uncond_pp_coverage_threespace(threespace):
#     cov = get_uncond_pp_coverage(threespace[0].trees)
#     assert cov == 1, "Coverage failed (3space)"


def test_get_uncond_pp_coverage_threespace(fourspace):
    cov = get_uncond_pp_coverage(fourspace[0].trees)
    cov == 1, "Coverage failed (4space)!"


def test_get_greedy_uncond_pp_tree(ten_taxa_multichain):
    dict_part = get_dict_of_unconditional_partitions(ten_taxa_multichain[0].trees)
    tree = get_greedy_uncond_pp_tree(dict_part, len(ten_taxa_multichain[0].trees[0]))
    assert True


# def test_get_tree_from_partition(ten_taxa_cMChain):
#     distance_sum = 0
#     for tree in ten_taxa_cMChain[0].trees:
#         cp = get_conditional_partitions(tree)
#         cp = [cp[k] for k in sorted(cp, reverse=True)[1:]]  # adapting for now
#         nt = get_tree_from_partition(cp, 10)
#         distance_sum += nt.fp_distance(tree)
#     assert distance_sum == 0, "Tree from partition failed!"
#
#
# def test_get_tree_from_partition20(twenty_taxa_tts):
#     cp = get_conditional_partitions(twenty_taxa_tts[0])
#     cp = [cp[k] for k in sorted(cp, reverse=True)[1:]]  # adapting for now
#     nt = get_tree_from_partition(cp, len(twenty_taxa_tts[0]))
#     assert nt.fp_distance(twenty_taxa_tts[0]) == 0, "Get tree form partition for 20 taxa failed!"
#
#
# def test_sample_from_dict_partition(ten_taxa_cMChain):
#     dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
#     t = sample_from_dict_partition(dict_part, samples=100, n_taxa=len(ten_taxa_cMChain[0].trees[0]))
#     assert True
#
#
# def test_search_maxpp_tree(ten_taxa_cMChain):
#     dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
#     t, t_lp = search_maxpp_tree(dict_part, n_taxa=len(ten_taxa_cMChain[0].trees[0]))
#     baseline = get_greedy_pp_tree(dict_part, len(ten_taxa_cMChain[0].trees[0]))
#     baseline_p = get_pp(baseline, dict_part, log=True)
#     assert t_lp > baseline_p, "Branch and bound search went wrong."
