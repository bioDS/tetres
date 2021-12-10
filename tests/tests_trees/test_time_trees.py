import ete3
import pytest
import warnings

from treeoclock.trees.time_trees import TimeTree, findpath_distance, findpath_path, get_mapping_dict,\
    TimeTreeSet, get_rank_neighbours, get_nni_neighbours, neighbourhood, nwk_to_cluster, DifferentNbrTaxa, free_tree_list
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


def test_findpath_path_typeerror():
    with pytest.raises(TypeError):
        findpath_path(1, 2)


def test_fp_path_warning_ctree(five_taxa_tts):
    t = five_taxa_tts
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        findpath_path(t[0].ctree, t[1].ctree)
        assert issubclass(w[-1].category, UserWarning)


def test_fp_path_warning_timetree(five_taxa_tts):
    t = five_taxa_tts
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        findpath_path(t[0], t[1])
        assert issubclass(w[-1].category, UserWarning)


def test_fp_path_warning_ete3(five_taxa_tts):
    t = five_taxa_tts
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        findpath_path(t[0].etree, t[1].etree)
        assert issubclass(w[-1].category, UserWarning)


def test_timetree_fp_path(five_taxa_newick_list, five_taxa_list_distances):
    out = []
    t = [TimeTree(i) for i in five_taxa_newick_list]
    for i in t:
        ext_path = []
        for j in t:
            cur_path = i.fp_path(j)
            ext_path.append(cur_path.num_trees)
            free_tree_list(cur_path)
        out.extend(ext_path)
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
            cur_path = findpath_path(i, j)
            out.append(cur_path.num_trees)
            free_tree_list(cur_path)
    assert out == five_taxa_list_distances, f"findpath_distances for TimeTree class wrong {out}"


def test_findpath_path_ctree(five_taxa_newick_list, five_taxa_list_distances):
    t = [ete3.Tree(i) for i in five_taxa_newick_list]
    out = []
    for i in t:
        for j in t:
            cur_path = findpath_path(i, j)
            out.append(cur_path.num_trees)
            free_tree_list(cur_path)
    assert out == five_taxa_list_distances, f"findpath_distances for ete3.Tree class wrong {out}"


def test_findpath_path_ete3tree(five_taxa_newick_list, five_taxa_list_distances):
    t = [ete3_to_ctree(ete3.Tree(i)) for i in five_taxa_newick_list]
    out = []
    for i in t:
        for j in t:
            cur_path = findpath_path(i, j)
            out.append(cur_path.num_trees)
            free_tree_list(cur_path)
    assert out == five_taxa_list_distances, f"findpath_distances for c-TREE class wrong {out}"


def test_get_mapping_dict(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    d = get_mapping_dict(f"{dir.path}/test.nex")
    f = open(f"{dir.path}/test.nex", "r")
    assert d == {1: 't1', 2: 't5', 3: "t4", 4: "t3", 5: "t2"}


def test_timetreeset_construction(dir, five_taxa_nexus_string):
    dir.write("test.nex", five_taxa_nexus_string)
    assert TimeTreeSet(f'{dir.path}/test.nex'), 'Construction not possible!'


def test_timetreeset_construction_filenotexists(dir):
    with pytest.raises(Exception):
        TimeTreeSet(f'{dir.path}/test5.nex')


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


def test_timetreeset_fp_path(five_taxa_tts, five_taxa_list_distances):
    t = five_taxa_tts
    out = []
    for i in range(len(t)):
        ext_list = []
        for j in range(len(t)):
            cur_path = t[i].fp_path(t[j])
            ext_list.append(cur_path.num_trees)
            free_tree_list(cur_path)
        out.extend(ext_list)
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
        five_taxa_tts.change_mapping({1: "ta1", 2: "ta2", 3: "ta3", 4: "ta4", 5: "ta5"})


def test_apply_new_taxa_map_difftaxa_more_error(five_taxa_tts):
    with pytest.raises(ValueError):
        five_taxa_tts[0].apply_new_taxa_map({1: "t1", 2: "t2", 3: "t3", 4: "t4", 5: "t5", 6: "t6"}, five_taxa_tts.map)


def test_apply_new_taxa_map_difftaxa_less_error(five_taxa_tts):
    with pytest.raises(ValueError):
        five_taxa_tts[0].apply_new_taxa_map({1: "t1", 4: "t2", 3: "t3", 5: "t5"}, five_taxa_tts.map)


def test_apply_new_taxa_map_difftaxanames_error(five_taxa_tts):
    with pytest.raises(ValueError):
        five_taxa_tts[0].apply_new_taxa_map({1: "ta1", 2: "ta2", 3: "ta3", 4: "ta4", 5: "ta5"}, five_taxa_tts.map)


def test_apply_new_taxa_map_same_map(five_taxa_tts):
    new_t = five_taxa_tts[0].apply_new_taxa_map(five_taxa_tts.map, five_taxa_tts.map)
    assert [new_t.fp_distance(five_taxa_tts[0]),
            new_t.fp_distance(five_taxa_tts[1]),
            new_t.fp_distance(five_taxa_tts[2])] == [0, 1, 2], 'Apply new map does not work!'


def test_apply_new_taxa_map(five_taxa_tts):
    new_t = five_taxa_tts[0].apply_new_taxa_map(old_map=five_taxa_tts.map,
                                                new_map={1: "t1", 2: "t2", 3: "t3", 4: "t4", 5: "t5"})
    assert new_t.fp_distance(five_taxa_tts[0]) == 5, "Apply new map failed!"


def test_change_mapping(five_taxa_tts, five_taxa_list_distances):
    new_trees = five_taxa_tts
    new_trees.change_mapping({1: "t1", 2: "t2", 3: "t3", 4: "t4", 5: "t5"})
    out = []
    for i in new_trees:
        out.extend([i.fp_distance(j) for j in new_trees])
    assert out == five_taxa_list_distances, f'Change_mapping failed!'


def test_change_mapping_more(twelve_taxa_tts, twelve_taxa_tts_start):
    new_start = twelve_taxa_tts_start
    new_start.change_mapping(twelve_taxa_tts.map)
    assert new_start[0].fp_distance(twelve_taxa_tts[2279]) == 0, "Change mapping failed!"

