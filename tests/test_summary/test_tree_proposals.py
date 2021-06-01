
from treeoclock.summary._tree_proposals import search_neighbourhood_greedy
from treeoclock.trees.time_trees import TimeTreeSet


def test_search_neighbourhood_greedy_ncorenone(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    res = search_neighbourhood_greedy(t=trees[1], trees=trees, t_value=10)
    assert trees[0].fp_distance(res) == 0, "Greedy neighbourhood search with n_cores=None failed!"


def test_search_neighbourhood_greedy_ncore1(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    res = search_neighbourhood_greedy(t=trees[1], trees=trees, t_value=10, n_cores=1)
    assert trees[0].fp_distance(res) == 0, "Greedy neighbourhood search with n_cores=1 failed!"


def test_search_neighbourhood_greedy_ncore4(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    res = search_neighbourhood_greedy(t=trees[1], trees=trees, t_value=10, n_cores=4)
    assert trees[0].fp_distance(res) == 0, "Greedy neighbourhood search with n_cores=4 failed!"
