import pytest

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
    assert sos == 288830, "Wrong SoS value for twelve taxa dataset!"


# Tests for inc_sub
def test_centroid_inc_sub(five_taxa_tts):
    cen, sos = Centroid(variation="inc_sub").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Greedy centroid variation failed!"


def test_centroid_inc_sub_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="inc_sub", start=10, subsample_size=2500).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Wrong SoS value for twelve taxa dataset!"


# Tests for iter_sub
def test_centroid_iter_sub(five_taxa_tts):
    cen, sos = Centroid(variation="iter_sub").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Greedy centroid variation failed!"


def test_centroid_iter_sub_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="iter_sub", start=10, subsample_size=2500).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Wrong SoS value for twelve taxa dataset!"


# Tests for separate
def test_centroid_separate(five_taxa_tts):
    cen, sos = Centroid(variation="separate").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Greedy centroid variation failed!"


def test_centroid_separate_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="separate", start=10).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Wrong SoS value for twelve taxa dataset!"


# Tests for onlyone
def test_centroid_onlyone(five_taxa_tts):
    cen, sos = Centroid(variation="onlyone").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Greedy centroid variation failed!"


def test_centroid_onlyone_12data(twelve_taxa_tts):
    # Fixing the subsample size to the number of trees in the set for testing
    cen, sos = Centroid(variation="onlyone", start=10).compute_centroid(twelve_taxa_tts)
    assert sos == 288830, "Wrong SoS value for twelve taxa dataset!"


def test_centroid_start_TTS_lengthError(twelve_taxa_tts):
    with pytest.raises(ValueError):
        Centroid(variation="inc_sub", n_cores=1, start=twelve_taxa_tts).compute_centroid(twelve_taxa_tts)


def test_centroid_start_TTS_taxaerror(twelve_taxa_tts, twenty_taxa_tts_start):
    with pytest.raises(ValueError):
        Centroid(variation="inc_sub", n_cores=1, start=twenty_taxa_tts_start).compute_centroid(twelve_taxa_tts)

# todo need test for when mapping is wrong
# test the general case if works
# test the error if wrong number of taxa
# test the error only one tree in the starting tree set

