from treeoclock.trees.time_trees import TimeTree

def test_fp_distance(five_taxa_newick_list):
    out = []
    t0 = TimeTree(five_taxa_newick_list[0])
    t1 = TimeTree(five_taxa_newick_list[1])
    out.append(t0.fp_distance(t0))
    out.append(t1.fp_distance(t1))
    out.append(t0.fp_distance(t1))
    out.append(t1.fp_distance(t0))
    assert out == [0, 0, 1, 1], f'fp_distance wrong {out} instead of [0,0,1,1]'


def test_fp_path():
    assert True


def test_get_newick():
    assert True
