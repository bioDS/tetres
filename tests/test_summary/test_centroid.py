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


def test_centroid_greedy(five_taxa_tts):
    cen, sos = Centroid(variation="greedy").compute_centroid(five_taxa_tts)
    t_cen = TimeTree('((3:2,(4:1,5:1):1):2,(1:3,2:3):1);')
    assert (cen.fp_distance(t_cen), sos) == (0, 5), "Greedy centroid variation failed!"


def test_centroid_compute_centroid_var_error2():
    with pytest.raises(ValueError):
        Centroid(variation="Greedy", n_cores=1).compute_centroid(TimeTreeSet())


def test_centroid_greedy_20data():
    TimeTreeSet(f"{Path(__file__).parent.absolute()}/data/20Taxa.trees")
    assert True
# TODO new variations testing, include a bigger dataset via a trees file
