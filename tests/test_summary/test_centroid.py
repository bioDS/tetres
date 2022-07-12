import pytest
import warnings
import random
import os

from pathlib import Path

from treeoclock.summary.centroid import Centroid
from treeoclock.trees.time_trees import TimeTreeSet, TimeTree


def test_centroid_construction_empty():
    assert Centroid(), "Empty construction failed!"


def test_centroid_construction():
    assert Centroid(variation="greedy", n_cores=4, select='last', start='first'), "Construction not possible!"


def test_centroid_compute_centroid_var_error1():
    with pytest.raises(ValueError):
        Centroid(variation="TEST", n_cores=1).compute_centroid(TimeTreeSet())


def test_centroid_compute_centroid_n_cores_error():
    with pytest.raises(ValueError):
        Centroid(variation="greedy", n_cores=0.9).compute_centroid(TimeTreeSet())


def test_centroid_compute_centroid_select_error():
    with pytest.raises(ValueError):
        Centroid(variation="TEST", n_cores=1, select='test').compute_centroid(TimeTreeSet())


def test_centroid_compute_centroid_start_error():
    with pytest.raises(ValueError):
        Centroid(variation="TEST", n_cores=1, start='test').compute_centroid(TimeTreeSet())


def test_centroid_compute_centroid_start_oor_error():
    with pytest.raises(ValueError):
        Centroid(variation="TEST", n_cores=1, start=1).compute_centroid(TimeTreeSet())


def test_centroid_compute_centroid_var_error2():
    with pytest.raises(ValueError):
        Centroid(variation="Greedy", n_cores=1).compute_centroid(TimeTreeSet())


# Tests for greedy
def test_centroid_greedy(five_taxa_tts):
    cen, sos = Centroid(variation="greedy").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Greedy centroid variation failed!"


def test_centroid_greedy_12data(twelve_taxa_tts):
    cen, sos = Centroid(variation="greedy", start=10).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Greedy: Wrong SoS value for twelve taxa dataset!"


# Tests for greedy_omp
def test_centroid_greedy(five_taxa_tts):
    cen, sos = Centroid(variation="greedy_omp").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Greedy_omp centroid variation failed!"


def test_centroid_greedy_12data(twelve_taxa_tts):
    cen, sos = Centroid(variation="greedy_omp", start=10).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Greedy_omp: Wrong SoS value for twelve taxa dataset!"


# Tests for inc_sub
def test_centroid_inc_sub(five_taxa_tts):
    cen, sos = Centroid(variation="inc_sub").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "IncSub centroid variation failed!"


def test_centroid_inc_sub_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="inc_sub", start=10, subsample_size=2500).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "IncSub: Wrong SoS value for twelve taxa dataset!"


def test_centroid_inc_sub_maxiter(five_taxa_tts):
    cen, sos = Centroid(variation="inc_sub", max_iterations=1).compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "IncSub centroid variation failed with max_iteration parameter!"


def test_centroid_inc_sub_12data_maxiter(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="inc_sub", start=10, subsample_size=2500, max_iterations=1).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "IncSub: Wrong SoS value for twelve taxa dataset!"


def test_centroid_inc_sub_12data_maxiter(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="inc_sub", start=10, subsample_size=2500, max_iterations=0).compute_centroid(twelve_taxa_tts)
    assert cen.fp_distance(twelve_taxa_tts[10]) == 0, "IncSub with max_iterations=0 failed!"


# Tests for iter_sub
def test_centroid_iter_sub(five_taxa_tts):
    cen, sos = Centroid(variation="iter_sub", start="sortFM").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "IterSub centroid variation failed!"


def test_centroid_iter_sub_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="iter_sub", start=10, subsample_size=2500).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "IterSub: Wrong SoS value for twelve taxa dataset!"


def test_centroid_iter_sub_12data_maxiter(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="iter_sub", start=10, subsample_size=2500, max_iterations=1).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "IterSub: Wrong SoS value for twelve taxa dataset!"


def test_centroid_iter_sub_12data_maxiter(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="iter_sub", start=10, subsample_size=2500, max_iterations=0).compute_centroid(twelve_taxa_tts)
    assert cen.fp_distance(twelve_taxa_tts[10]) == 0, "IncSub with max_iterations=0 failed!"


# Tests for separate
def test_centroid_separate(five_taxa_tts):
    cen, sos = Centroid(variation="separate").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Separate centroid variation failed!"


def test_centroid_separate_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="separate", start=10).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Separate: Wrong SoS value for twelve taxa dataset!"


# Tests for onlyone
def test_centroid_onlyone(five_taxa_tts):
    cen, sos = Centroid(variation="onlyone").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Onlyone centroid variation failed!"


def test_centroid_onlyone_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="onlyone", start=10).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Onlyone: Wrong SoS value for twelve taxa dataset!"


# Tests for update_with_one
def test_centroid_updateone(five_taxa_tts):
    cen, sos = Centroid(variation="update_with_one").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "UpdateOne centroid variation failed!"


def test_centroid_updateone_12data(twelve_taxa_tts):
    state = random.getstate()  # get the random seed state
    random.seed(28)  # Fixing the seed to get the same result
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="update_with_one", start=10).compute_centroid(twelve_taxa_tts)
    random.setstate(state)  # reset the random seed to previous state
    assert sos == 288830, "UpdateOne: Wrong SoS value for twelve taxa dataset!"


# Tests for online
def test_centroid_online(five_taxa_tts):
    cen, sos = Centroid(variation="online").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Online centroid variation failed!"


def test_centroid_online_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="online", start=10).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Online: Wrong SoS value for twelve taxa dataset!"


# more tests
def test_centroid_start_tts_wrong_taxa_error(twelve_taxa_tts, twenty_taxa_tts_start):
    with pytest.raises(ValueError):
        Centroid(variation="inc_sub", n_cores=1, start=twenty_taxa_tts_start).compute_centroid(twelve_taxa_tts)


def test_centroid_start_tts_wrong_map_check_map(twelve_taxa_tts, twelve_taxa_tts_start):
    my_cen = Centroid(variation="inc_sub", n_cores=1, start=twelve_taxa_tts_start)
    assert my_cen.start[0].fp_distance(twelve_taxa_tts[2279]) == 0,\
        "Centroid with TTS start different map failed to change map!"


def test_centroid_start_tts_wrong_map_running(twelve_taxa_tts, twelve_taxa_tts_start):
    my_cen = Centroid(variation="inc_sub", n_cores=1, start=twelve_taxa_tts_start)
    assert my_cen.compute_centroid(twelve_taxa_tts),\
        "Centroid start with TTS different map failed!"


def test_centroid_treelogfile_wrong_path(twelve_taxa_tts):
    my_cen = Centroid(tree_log_file=f'{Path(__file__).parent.absolute()}/data/twelve_log.trees')
    with pytest.raises(FileNotFoundError):
        my_cen.compute_centroid(twelve_taxa_tts)


def test_centroid_treelogfile_warning(twelve_taxa_tts):
    my_cen = Centroid(tree_log_file=f'{Path(__file__).parent.parent.absolute()}/data/logfile_12.trees')
    if not os.path.exists(f'{Path(__file__).parent.parent.absolute()}/data/logfile_12.trees'):
        open(f'{Path(__file__).parent.parent.absolute()}/data/logfile_12.trees', "x")
    with pytest.warns(UserWarning):
        my_cen.compute_centroid(twelve_taxa_tts)


# todo test for it working
def test_centroid_treelogfile_greedy(twelve_taxa_tts):
    my_cen = Centroid(tree_log_file=f'{Path(__file__).parent.parent.absolute()}/data/logfile_12.trees')
    with warnings.catch_warnings():
        # Catching the warning so that it is not displayed for testing
        warnings.filterwarnings("ignore", category=UserWarning)
        assert my_cen.compute_centroid(twelve_taxa_tts)


def test_centroid_treelogfile_incsub(twelve_taxa_tts):
    my_cen = Centroid(variation="inc_sub",
                      tree_log_file=f'{Path(__file__).parent.parent.absolute()}/data/logfile_12.trees')
    with warnings.catch_warnings():
        # Catching the warning so that it is not displayed for testing
        warnings.filterwarnings("ignore", category=UserWarning)
        assert my_cen.compute_centroid(twelve_taxa_tts)


