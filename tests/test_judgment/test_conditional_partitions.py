from treeoclock.judgment.conditional_partitions import get_conditional_partitions, get_dict_of_partitions, get_pp, get_pp_coverage


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
    assert True


def test_get_pp_coverage_threespace(fourspace):
    cov = get_pp_coverage(fourspace[0].trees)
    assert True
