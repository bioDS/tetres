import pytest

from treeoclock.summary.centroid import Centroid
from treeoclock.trees.time_trees import TimeTreeSet


def test_centroid_construction_empty():
    assert Centroid(), "Empty construction failed!"


def test_centroid_construction():
    assert Centroid(variation="inc-sub", n_cores=None, start='first'), "Construction not possible!"


def test_centroid_compute_centroid_var_error():
    with pytest.raises(ValueError):
        Centroid(variation="TEST", n_cores=1)


def test_centroid_compute_centroid_n_cores_error():
    with pytest.raises(ValueError):
        Centroid(variation="greedy", n_cores=0.6)


def test_centroid_greedy():
    # TODO Needs to be adapted
    assert Centroid(variation="greedy", n_cores=4).compute_centroid(TimeTreeSet()) == 0, "Greedy Variation Failed!"


def test_centroid_n_cores_none():
    # TODO Needs to be adapted
    assert Centroid(variation="greedy").compute_centroid(TimeTreeSet()) == 0, "Greedy with n_cores=None failed!"


