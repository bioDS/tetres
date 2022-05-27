import pytest
from pathlib import Path
from treeoclock.judgment.mchain import coupled_MChains
import os
from treeoclock.visualize.graphs import *


def test_visualize_matrix(ten_taxa_cMChain):
    visualize_largest_cc(ten_taxa_cMChain, 0)
    assert True



