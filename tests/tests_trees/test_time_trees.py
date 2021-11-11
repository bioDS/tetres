import ete3
import pytest

from treeoclock.trees.time_trees import TimeTree, findpath_distance, findpath_path, get_mapping_dict,\
    TimeTreeSet, get_rank_neighbours, get_nni_neighbours, neighbourhood, nwk_to_cluster, DifferentNbrTaxa
from treeoclock.trees._converter import ete3_to_ctree


def test_timetree_construction(five_taxa_newick_list):
    assert TimeTree(five_taxa_newick_list[0]), "Construction not possible!"


def test_timetree_fp_distance(five_taxa_newick_list, five_taxa_list_distances):
    out = []
    t = [TimeTree(i) for i in five_taxa_newick_list]
    for i in t:
        out.extend([i.fp_distance(j) for j in t])
    assert out == five_taxa_list_distances, f'fp_distance wrong {out}'


def test_timetree_fp_distance_ctree_differentnbrtaxa(five_taxa_tts, twelve_taxa_tts):
    with pytest.raises(DifferentNbrTaxa):
        findpath_distance(five_taxa_tts[0].ctree, twelve_taxa_tts[0].ctree)


def test_timetree_fp_distance_etree_differentnbrtaxa(five_taxa_tts, twelve_taxa_tts):
    with pytest.raises(DifferentNbrTaxa):
        findpath_distance(five_taxa_tts[0].etree, twelve_taxa_tts[0].etree)


def test_timetree_fp_distance_timetree_differentnbrtaxa(five_taxa_tts, twelve_taxa_tts):
    with pytest.raises(DifferentNbrTaxa):
        findpath_distance(five_taxa_tts[0], twelve_taxa_tts[0])


# TODO catch the error for the fp_path
#  also add the test for all different data types that could be used


def test_timetree_fp_path(five_taxa_newick_list, five_taxa_list_distances):
    out = []
    t = [TimeTree(i) for i in five_taxa_newick_list]
    for i in t:
        out.extend([len(i.fp_path(j)) for j in t])
    assert out == five_taxa_list_distances, f"fp_path distances wrong {out}"


def test_timetree_get_newick(five_taxa_newick_list):
    t = [TimeTree(i) for i in five_taxa_newick_list]
    out = []
    for i in t:
        out.append(i.get_newick())
    assert out == five_taxa_newick_list, "TimeTree does not return correct newick strings"


def test_timetree_copy(five_taxa_newick_list):
    t = [TimeTree(i) for i in five_taxa_newick_list]
    out = [i.copy() for i in t]
    bout = []
    for i in range(len(out)):
        bout.append(out[i].get_newick() == t[i].get_newick())
    assert bout.count(False) == 0, f"Copying a TimeTree failed!"


def test_timetree_neighbours(five_taxa_newick_list, five_taxa_0_all_neighbours):
    t = TimeTree(five_taxa_newick_list[0])
    n = t.neighbours()
    assert set([i.get_newick() for i in n]) == set(five_taxa_0_all_neighbours), "TimeTree neighbours wrong!"


def test_timetree_rank_neighbours(five_taxa_newick_list):
    t = TimeTree(five_taxa_newick_list[0])
    n = t.rank_neighbours()
    assert n[0].get_newick() == "((4:3,(3:1,2:1):2):1,(1:2,5:2):2);", "TimeTree rank neighbour wrong!"

def test_timetree_nni_neighbours(five_taxa_newick_list, five_taxa_0_nni_neighbours):
    t = TimeTree(five_taxa_newick_list[0])
    n = t.nni_neighbours()
    assert set([i.get_newick() for i in n]) == set(five_taxa_0_nni_neighbours), "TimeTree nni neighbours wrong!"


def test_findpath_distance_timetree(five_taxa_newick_list, five_taxa_list_distances):
    t = [TimeTree(i) for i in five_taxa_newick_list]
    out = []
    for i in t:
        for j in t:
            out.append(findpath_distance(i, j))
    assert out == five_taxa_list_distances, f"findpath_distances for TimeTree class wrong {out}"


def test_findpath_distance_ctree(five_taxa_newick_list, five_taxa_list_distances):
    t = [ete3.Tree(i) for i in five_taxa_newick_list]
    out = []
    for i in t:
        for j in t:
            out.append(findpath_distance(i, j))
    assert out == five_taxa_list_distances, f"findpath_distances for ete3.Tree class wrong {out}"


def test_findpath_distance_ete3tree(five_taxa_newick_list, five_taxa_list_distances):
    t = [ete3_to_ctree(ete3.Tree(i)) for i in five_taxa_newick_list]
    out = []
    for i in t:
        for j in t:
            out.append(findpath_distance(i, j))
    assert out == five_taxa_list_distances, f"findpath_distances for c-TREE class wrong {out}"


def test_findpath_path_timetree(five_taxa_newick_list, five_taxa_list_distances):
    t = [TimeTree(i) for i in five_taxa_newick_list]
    out = []
    for i in t:
        for j in t:
            out.append(len(findpath_path(i, j)))
    assert out == five_taxa_list_distances, f"findpath_distances for TimeTree class wrong {out}"


def test_findpath_path_ctree(five_taxa_newick_list, five_taxa_list_distances):
    t = [ete3.Tree(i) for i in five_taxa_newick_list]
    out = []
    for i in t:
        for j in t:
            out.append(len(findpath_path(i, j)))
    assert out == five_taxa_list_distances, f"findpath_distances for ete3.Tree class wrong {out}"


def test_findpath_path_ete3tree(five_taxa_newick_list, five_taxa_list_distances):
    t = [ete3_to_ctree(ete3.Tree(i)) for i in five_taxa_newick_list]
    out = []
    for i in t:
        for j in t:
            out.append(len(findpath_path(i, j)))
    assert out == five_taxa_list_distances, f"findpath_distances for c-TREE class wrong {out}"


def test_get_mapping_dict(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    d = get_mapping_dict(f"{dir.path}/test.nex")
    f = open(f"{dir.path}/test.nex", "r")
    assert d == {1: 't1', 2: 't5', 3: "t4", 4: "t3", 5: "t2"}


def test_timetreeset_construction(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    assert TimeTreeSet(f'{dir.path}/test.nex'), 'Construction not possible!'


def test_timetreeset_construction_empty():
    assert type(TimeTreeSet()) is TimeTreeSet, "Empty construction failed!"


def test_timetreeset_len(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    t = TimeTreeSet(f'{dir.path}/test.nex')
    assert len(t) == 3, f'The nexus file should contain 3 trees, but only {len(t)} are being returned!'


def test_timetreeset_fp_distance(dir, five_taxa_nexus_string, five_taxa_list_distances):
    dir.write("test.nex", five_taxa_nexus_string)
    t = TimeTreeSet(f'{dir.path}/test.nex')
    out = []
    for i in range(len(t)):
        out.extend([t[i].fp_distance(t[j]) for j in range(len(t))])
    assert out == five_taxa_list_distances, 'TimeTreeSet fp_distance() wrong!'


def test_timetreeset_fp_path(dir, five_taxa_nexus_string, five_taxa_list_distances):
    dir.write("test.nex", five_taxa_nexus_string)
    t = TimeTreeSet(f'{dir.path}/test.nex')
    out = []
    for i in range(len(t)):
        out.extend([len(t[i].fp_path(t[j])) for j in range(len(t))])
    assert out == five_taxa_list_distances, 'TimeTreeSet fp_path() wrong!'


def test_get_rank_neighbours(five_taxa_newick_list):
    t = TimeTree(five_taxa_newick_list[0])
    rn = get_rank_neighbours(t)
    assert rn[0].get_newick() == "((4:3,(3:1,2:1):2):1,(1:2,5:2):2);", "Rank neighbour wrong!"


def test_get_rank_neighbours2(seventeen_taxa_tree_newick, seventeen_taxa_rank_neighbours_newick):
    t = TimeTree(seventeen_taxa_tree_newick)
    rn = get_rank_neighbours(t)
    assert set(seventeen_taxa_rank_neighbours_newick) == set([i.get_newick(f=5) for i in rn]), "Rank neighbours wrong!"


def test_get_nni_neighbours(five_taxa_newick_list, five_taxa_0_nni_neighbours):
    t = TimeTree(five_taxa_newick_list[0])
    nn = get_nni_neighbours(t)
    assert set([i.get_newick() for i in nn]) == set(five_taxa_0_nni_neighbours), "NNI neighbours wrong!"


def test_get_nni_neighbours2(seventeen_taxa_tree_newick, seventeen_taxa_nni_neighbours_newick):
    t = TimeTree(seventeen_taxa_tree_newick)
    nn = get_nni_neighbours(t)
    assert set(seventeen_taxa_nni_neighbours_newick) == set([i.get_newick(f=5) for i in nn]), "NNI neighbours wrong!"


def test_neighbourhood(five_taxa_newick_list, five_taxa_0_all_neighbours):
    t = TimeTree(five_taxa_newick_list[0])
    n = neighbourhood(t)
    assert set([i.get_newick() for i in n]) == set(five_taxa_0_all_neighbours), "Neighbourhood wrong!"


def test_neighbourhood2(seventeen_taxa_tree_newick, seventeen_taxa_nni_neighbours_newick,
                        seventeen_taxa_rank_neighbours_newick):
    t = TimeTree(seventeen_taxa_tree_newick)
    n = neighbourhood(t)
    assert set([i.get_newick() for i in n]) == \
           set(seventeen_taxa_nni_neighbours_newick+seventeen_taxa_rank_neighbours_newick), "Neighbourhood wrong!"


def test_get_clades(five_taxa_newick_list):
    t = TimeTree(five_taxa_newick_list[0])
    out = t.get_clades()
    assert out == {frozenset({"1", "5"}), frozenset({"2", "3"}), frozenset({"2", "3", "4"})}, "Function get_clades failed!"


def test_nwk_to_cluster():
    assert nwk_to_cluster("((A:1,B:1):1, C:2);"), "Function nwk_to_cluster failed!"


def test_nwk_to_cluster_exception1():
    # No semicolon in the end
    with pytest.raises(Exception):
        nwk_to_cluster("((A:1,B:1):1, C:2)")


def test_nwk_to_cluster_exception2():
    # No bracket at position -2
    with pytest.raises(Exception):
        nwk_to_cluster("((A:1,B:1):1, C:2;")


def test_nwk_to_cluster_exception3():
    # To many opening brackets (
    with pytest.raises(Exception):
        nwk_to_cluster("(((A:1,B:1):1, C:2);")


def test_nwk_to_cluster_exception4():
    # To many closing brackets )
    with pytest.raises(Exception):
        nwk_to_cluster("((A:1,B:1)):1, C:2);")


def test_get_common_clades(five_taxa_tts):
    clades = five_taxa_tts.get_common_clades()
    assert clades == set([frozenset(["2", "1"]), frozenset(["3", "4", "5"])])


def test_change_mapping_samemap_error(five_taxa_tts):
    with pytest.raises(ValueError):
        five_taxa_tts.change_mapping(five_taxa_tts.map)


def test_change_mapping_difftaxa_more_error(five_taxa_tts):
    with pytest.raises(ValueError):
        five_taxa_tts.change_mapping({1: "t1", 2: "t2", 3: "t3", 4: "t4", 5: "t5", 6: "t6"})


def test_change_mapping_difftaxa_less_error(five_taxa_tts):
    with pytest.raises(ValueError):
        five_taxa_tts.change_mapping({1: "t1", 4: "t2", 3: "t3", 5: "t5"})


def test_change_mapping_difftaxanames_error(five_taxa_tts):
    with pytest.raises(ValueError):
        five_taxa_tts.change_mapping({1: "t1", 2: "t2", 3: "t3", 4: "t4", 5: "t5"})


def test_change_mapping(five_taxa_tts):
    # todo!!!
    assert five_taxa_tts.change_mapping({1: "t1", 2: "t2", 3: "t3", 4: "t4", 5: "t5"})


def test_apply_new_taxa_map():
    # Todo
    assert True, "apply_new_taxa_map failed!"


def test_change_mapping():
    # Todo
    assert True, "change_mapping failed!"
