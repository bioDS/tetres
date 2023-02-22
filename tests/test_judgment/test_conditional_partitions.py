from treeoclock.judgment.conditional_partitions import get_conditional_partitions, get_dict_of_partitions, get_pp, get_pp_coverage, get_greedy_pp_tree, get_tree_from_partition
from treeoclock.trees._converter import ctree_to_ete3

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


def test_get_pp_coverage(ten_taxa_cMChain):
    cov = get_pp_coverage(ten_taxa_cMChain[0].trees)
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


def test_get_tree_from_partition(ten_taxa_cMChain):
    distance_sum = 0
    for tree in ten_taxa_cMChain[0].trees:
        cp = get_conditional_partitions(tree)
        cp = [cp[k] for k in sorted(cp, reverse=True)[1:]]  # adapting for now
        nt = get_tree_from_partition(cp, 10)
        distance_sum += nt.fp_distance(tree)
    assert distance_sum == 0, "Tree from partition failed!"
