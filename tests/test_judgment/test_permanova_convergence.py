import pytest
from treeoclock.judgment.permanova_convergence import *
import os


def test_permanova_cut(ten_taxa_cMChain):
    s, e = permanova_cut(ten_taxa_cMChain, 0, 1)
    assert (s, e) == (1, 10), "Permanove Failed!"


def test_permanova_pvalue_plot(ten_taxa_cMChain):
    permanova_pvalue_plot(ten_taxa_cMChain, 0, 1, smoothing_l=[0.9, 0.5])
    assert True
