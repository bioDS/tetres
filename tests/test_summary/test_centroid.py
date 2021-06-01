import pytest

from treeoclock.summary.centroid import Centroid
from treeoclock.trees.time_trees import TimeTreeSet


def test_centroid_construction_empty():
    assert Centroid(), "Empty construction failed!"


def test_centroid_construction():
    assert Centroid(variation="inc-sub"), "Construction not possible!"


def test_centroid_compute_centroid_var_error():
    with pytest.raises(ValueError):
        Centroid(variation="TEST", n_cores=1).compute_centroid(TimeTreeSet())


def test_centroid_compute_centroid_n_cores_error():
    with pytest.raises(ValueError):
        Centroid(variation="greedy", n_cores=0.9).compute_centroid(TimeTreeSet())


def test_centroid_greedy():
    # TODO Needs to be adapted
    assert Centroid(variation="greedy").compute_centroid(TimeTreeSet()) == 0, "Greedy Variation Failed"




