
from treeoclock.summary.compute_sos import compute_sos
from treeoclock.trees.time_trees import TimeTree, TimeTreeSet


def test_compute_sos_one_core(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos(t=trees[0], trees=trees, n_cores=1) == 5, "Compute_sos failed, should be 0!"


def test_compute_sos_none_core(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos(t=trees[0], trees=trees, n_cores=None) == 5, "Compute_sos failed, should be 0!"


def test_compute_sos_four_cores(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos(t=trees[0], trees=trees, n_cores=4) == 5, "Compute_sos failed, should be 0!"