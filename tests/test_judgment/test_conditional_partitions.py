from tetres.judgment.conditional_partitions import get_conditional_partitions, get_dict_of_partitions, get_pp, get_greedy_pp_tree, get_tree_from_partition, sample_from_dict_partition, search_maxpp_tree, get_greedy_relaxed_pp_tree, are_compatible, convert_to_relaxed_partition_dict


def test_get_conditional_partitions(ten_taxa_cMChain):
    get_conditional_partitions(ten_taxa_cMChain[0].trees[0])
    assert True


def test_get_list_of_partitions(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    assert True


def test_get_pp(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    pp = get_pp(ten_taxa_cMChain[0].trees[0], dict_part)
    assert True


def test_get_pp_coverage_threespace(threespace):
    cov = get_pp_coverage(threespace[0].trees)
    assert cov == 1, "Coverage failed (3space)"


def test_get_pp_coverage_threespace(fourspace):
    cov = get_pp_coverage(fourspace[0].trees)
    cov == 1, "Coverage failed (4space)!"


def test_get_greedy_pp_tree(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    tree = get_greedy_pp_tree(dict_part, len(ten_taxa_cMChain[0].trees[0]))
    assert True


def test_are_compatible(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    keys = list(dict_part.keys())
    assert are_compatible(keys[1], keys[2]), "Are compatible function failed"


def test_are_compatible2(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    keys = list(dict_part.keys())
    assert not are_compatible(keys[1], keys[3]), "Are compatible function failed"


def test_get_tree_from_partition(ten_taxa_cMChain):
    distance_sum = 0
    for tree in ten_taxa_cMChain[0].trees:
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


def test_sample_from_dict_partition(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    t = sample_from_dict_partition(dict_part, samples=100, n_taxa=len(ten_taxa_cMChain[0].trees[0]))
    assert True


def test_search_maxpp_tree(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    t, t_lp = search_maxpp_tree(dict_part, n_taxa=len(ten_taxa_cMChain[0].trees[0]))
    baseline = get_greedy_pp_tree(dict_part, len(ten_taxa_cMChain[0].trees[0]))
    baseline_p = get_pp(baseline, dict_part, log=True)
    assert t_lp > baseline_p, "Branch and bound search went wrong."


def test_get_greedy_relaxed_pp_tree(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    relaxed_tree = get_greedy_relaxed_pp_tree(dict_part, len(ten_taxa_cMChain[0].trees[0]))
    greedy_tree = get_greedy_pp_tree(dict_part, len(ten_taxa_cMChain[0].trees[0]))
    assert relaxed_tree.fp_distance(greedy_tree) == 1, "Something went wrong with relaxed greedy!"


def test_convert_to_relaxed_partition_dict(ten_taxa_cMChain):
    dict_part = get_dict_of_partitions(ten_taxa_cMChain[0].trees)
    relaxed_dict = convert_to_relaxed_partition_dict(dict_part)
    greedy_tree = get_greedy_pp_tree(dict_part, len(ten_taxa_cMChain[0].trees[0]))
    relaxed_greedy_tree = get_greedy_pp_tree(relaxed_dict, len(ten_taxa_cMChain[0].trees[0]))
    relaxed_tree = get_greedy_relaxed_pp_tree(dict_part, len(ten_taxa_cMChain[0].trees[0]))
    assert (greedy_tree.fp_distance(relaxed_tree), greedy_tree.fp_distance(relaxed_greedy_tree), relaxed_greedy_tree.fp_distance(relaxed_tree)) == (1, 1, 0), "Converting dict of partitions to relaxed dict of partitions failed!"
