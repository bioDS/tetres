import pytest
from treeoclock.summary.centroid import Centroid
from treeoclock.summary.annotate_centroid import annotate_centroid
from treeoclock.trees.time_trees import TimeTree


def test_annotate_centroid(twelve_taxa_tts):
    cen, sos = Centroid(variation="greedy", start=10).compute_centroid(twelve_taxa_tts)
    acen = annotate_centroid(cen, twelve_taxa_tts)
    reversed_cen = TimeTree(acen.get_newick(f=5))
    assert acen.fp_distance(reversed_cen) == 0, "Centroid annotation failed"


# todo easy test case to see if it is doing the right thing with extracting heights of the ranks!