import pytest
from treeoclock.judgment.mchain import coupled_MChains
from treeoclock.judgment._gelman_rubin_diag import gelman_rubin_cut
from treeoclock.judgment._fast_gr_approximation import fast_grd_cut, gelman_rubin_fast_approx_parameter_choice_plot
import os


def test_original_values(ten_taxa_cMChain):
    # Testing the parameter choice plot function
    start, end = gelman_rubin_cut(ten_taxa_cMChain, 0, 1, smoothing=0.5, threshold_percentage=0.2,
                                  ess_threshold=0, pseudo_ess_range=100, smoothing_average="median")
    # ten_taxa_cMChain.gelman_rubin_parameter_choice_plot(i=0, j=1)
    assert [start, end] == [255, 318]


def test_fast_grd(ten_taxa_cMChain):
    start, end = fast_grd_cut(ten_taxa_cMChain, 0, 1, smoothing=0.5, threshold_percentage=0.2,
                                ess_threshold=0, pseudo_ess_range=100, smoothing_average="median")
    assert [start, end] == [-1, -1]

def test_cMChain_gelman_rubin_trace(ten_taxa_cMChain):
    # Testing the parameter choice plot function
    gelman_rubin_fast_approx_parameter_choice_plot(ten_taxa_cMChain, i=0, j=1)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_grd_fast_approximation_parameter_choices_ess-200_mean.png")
