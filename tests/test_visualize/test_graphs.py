import pytest
from pathlib import Path
from tetres.judgement.mchain import coupled_MChains
import os
from tetres.visualize.graphs import *


def test_visualize_matrix(ten_taxa_cMChain):
    visualize_largest_cc(ten_taxa_cMChain, 0)
    assert True



