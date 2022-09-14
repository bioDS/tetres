import pytest
from treeoclock.judgment.mchain import coupled_MChains
import os


def test_cMChain_gelman_rubin_plot(ten_taxa_cMChain):
    ten_taxa_cMChain.gelman_rubin_like_diagnostic_plot(samples="all")
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_grd_plot_all.png"), "Gelman Rubin Plot failed!"


def test_cMChain_gelman_rubin_trace(ten_taxa_cMChain):
    ten_taxa_cMChain.gelman_rubin_trace_plot(i=0, j=1)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_grd_singlevalue_evaluation.png"), \
        "Gelman Rubin Trace Plot failed!"


def test_cMChain_gress(ten_taxa_cMChain):
    ten_taxa_cMChain.gelman_rubin_trace_ess_plot(i=0, j=1, _overwrite=True)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_gress_arviz.png"), \
        "Gelman Rubin Trace Plot with ess failed!"


def test_cMChain_gress_with_ess_200(ten_taxa_cMChain):
    ten_taxa_cMChain.gelman_rubin_trace_ess_plot(i=0, j=1, ess=200, _overwrite=True)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_gress_200_arviz.png"), \
        "Gelman Rubin Trace Plot with ess failed!"


def test_cMChain_gress_with_ess_500(ten_taxa_cMChain):
    ten_taxa_cMChain.gelman_rubin_trace_ess_plot(i=0, j=1, ess=500, _overwrite=True)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_gress_500_arviz.png"), \
        "Gelman Rubin Trace Plot with ess failed!"
