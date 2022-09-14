import pytest
from treeoclock.judgment.mchain import coupled_MChains
import os


def test_cMChain_gelman_rubin_trace(ten_taxa_cMChain):
    # Testing the parameter choice plot function
    ten_taxa_cMChain.gelman_rubin_parameter_choice_plot(i=0, j=1)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_grd_singlevalue_evaluation.png"), \
        "Gelman Rubin Trace Plot failed!"


def test_cMChain_gelman_rubin_plot(ten_taxa_cMChain):
    # todo this test should delete the plot before recomputing it
    # todo test different sample parameters
    ten_taxa_cMChain.gelman_rubin_all_chains_density_plot(samples='all')
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_grd_density_full_chain_all.png"), "Gelman Rubin Plot failed!"


def test_cMChain_gress(ten_taxa_cMChain):
    # todo test different setting also
    ten_taxa_cMChain.gelman_rubin_cut(i=0, j=1, _overwrite=True)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/data/{ten_taxa_cMChain.name}_{0}_{1}_gelman_rubin_cutoff_no-esst_smoothing-0.5_thresholdp-0.5"), \
        "Gelman Rubin Trace Plot with ess failed!"
