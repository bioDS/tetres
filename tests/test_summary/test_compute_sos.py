
from treeoclock.summary.compute_sos import compute_sos, compute_sos_mt, compute_sos_omp
from treeoclock.trees.time_trees import TimeTreeSet


def test_compute_sos_one_core(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert [compute_sos_mt(t=ct, trees=trees, n_cores=1) for ct in trees] == [5,10,13], "Compute_sos_mp failed (1=cores), should be 5!"


def test_compute_sos_none_core(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos_mt(t=trees[0], trees=trees, n_cores=None) == 5, "Compute_sos_mp failed (None=cores), should be 5!"


def test_compute_sos_four_cores(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos_mt(t=trees[0], trees=trees, n_cores=4) == 5, "Compute_sos_mp failed (4=cores), should be 5!"

def test_compute_sos_omp_one_core(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos_omp(t=trees[0], trees=trees, n_cores=1) == 5, "Compute_sos_omp failed (1=cores), should be 5!"


def test_compute_sos_omp_none_core(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos_omp(t=trees[0], trees=trees, n_cores=None) == 5, "Compute_sos_omp failed (None=cores), should be 5!"


def test_compute_sos_omp_four_cores(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos_omp(t=trees[0], trees=trees, n_cores=4) == 5, "Compute_sos_omp failed (4=cores), should be 5!"


def test_compute_sos(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    trees = TimeTreeSet(f'{dir.path}/test.nex')
    assert compute_sos(trees[0], trees) == 5, f"Compute_sos failed, should be 5!"



