from treeoclock.trees.time_trees import TimeTree


def test_fp_distance(five_taxa_newick_list):
    out = []
    t = [TimeTree(i) for i in five_taxa_newick_list]
    for i in t:
        out.extend([i.fp_distance(j) for j in t])
    assert out == [0, 1, 2, 1, 0, 3, 2, 3, 0], f'fp_distance wrong {out}'


def test_fp_path():
    assert True


def test_get_newick():
    assert True
