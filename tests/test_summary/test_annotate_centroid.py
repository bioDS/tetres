import pytest
from treeoclock.summary.centroid import Centroid
from treeoclock.summary.annotate_centroid import annotate_centroid
from treeoclock.trees.time_trees import TimeTree


def test_annotate_centroid(twelve_taxa_tts):
    cen, sos = Centroid(variation="greedy", start=10).compute_centroid(twelve_taxa_tts)
    acen = annotate_centroid(cen, twelve_taxa_tts)
    reversed_cen = TimeTree(acen.get_newick(f=5))
    assert acen.fp_distance(reversed_cen) == 0, "Centroid annotation failed"


def test_annotate_tree(five_taxa_tts):
    t = annotate_centroid(five_taxa_tts[0], five_taxa_tts)
    assert t.get_newick(f=5) == '((1:3,2:3):1.66667,(3:2,(4:1,5:1):1):2.66667);',\
        "Annotation of tree with branch-lengths failed!"