import pytest
import os
from treeoclock.visualize.open_book import *
import matplotlib.pyplot as plt


def test_get_tree_coordinate(twenty_taxa_tts):
    t, h, bl = get_tree_coordinates(twenty_taxa_tts[0], [["1"], ["7"], ["9"]])
    assert [t, h, bl] == ["a|bc", 0.043759713940639045, 0.014579381202964016], "Get tree coordinates failed!"


def test_get_trees_embedding(twenty_taxa_tts):
    page_coords = get_trees_embedding(twenty_taxa_tts, [["1"], ["7"], ["9"]])
    assert len(page_coords["ab|c"]) == 8, "Something failed"


def test_plot_book(twenty_taxa_tts):
    # todo figure out a better triple example for this dataset
    fig = plot_book(twenty_taxa_tts, [["1"], ["7"], ["9"]])
    plt.show()


def test_get_all_triples(twenty_taxa_tts):
    # todo some test to check it works properly
    get_all_triples(twenty_taxa_tts[0].etree)


# def test_plot_pages(twenty_taxa_tts):
#     page_coords = get_trees_embedding(twenty_taxa_tts, [["1"], ["7"], ["9"]])
#     plot_pages(page1=page_coords["ab|c"], page2=page_coords ["bc|a"], page1_label="17|9", page2_label="79|1")
#     plot_pages(page1=page_coords["ab|c"], page2=page_coords["ac|b"], page1_label="17|9", page2_label="19|7")



