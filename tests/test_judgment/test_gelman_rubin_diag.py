import pytest
from tetres.judgment.mchain import coupled_MChains
import os


def test_cMChain_gelman_rubin_trace(ten_taxa_cMChain):
    # Testing the parameter choice plot function
    ten_taxa_cMChain.gelman_rubin_parameter_choice_plot(i=0, j=1)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_grd_parameter_choices_ess-list_mean.png"), \
        "Gelman Rubin Trace Plot failed!"


def test_cMChain_gelman_rubin_trace_subsampling(ten_taxa_cMChain):
    # Testing the parameter choice plot function
    ten_taxa_cMChain.gelman_rubin_parameter_choice_plot(i=0, j=1, _subsampling=1)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_grd_parameter_choices_subsampling-{1}_ess-list_mean.png"), \
        "Gelman Rubin Trace Plot failed!"


def test_cMChain_gelman_rubin_plot(ten_taxa_cMChain):
    # todo this test should delete the plot before recomputing it
    # todo test different sample parameters
    ten_taxa_cMChain.gelman_rubin_all_chains_density_plot(samples='all')
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_grd_density_full_chain_all.png"), "Gelman Rubin Plot failed!"


def test_cMChain_gr_cut(ten_taxa_cMChain):
    for _gr_boundary in [0.02, 0.01, 0.05]:
        for sm_avg in ['mean', 'median']:
            ten_taxa_cMChain.gelman_rubin_cut(i=0, j=1, _overwrite=True, _gr_boundary=_gr_boundary, smoothing_average=sm_avg)
            assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/data/{ten_taxa_cMChain.name}_{0}_{1}_gelman_rubin_cutoff_smoothing-0.5_{sm_avg}_boundary-{_gr_boundary}"), \
                f"Gelman Rubin Cut function failed! {_gr_boundary}-{sm_avg}"


def test_cMChain_gr_cut_subsample(ten_taxa_cMChain):
    for _gr_boundary in [0.02, 0.01, 0.05]:
        for sm_avg in ['mean', 'median']:
            ten_taxa_cMChain.gelman_rubin_cut(i=0, j=1, _overwrite=True, _subsampling=2, _gr_boundary=_gr_boundary, smoothing_average=sm_avg)
            assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/data/{ten_taxa_cMChain.name}_{0}_{1}_gelman_rubin_cutoff_subsampling-2_smoothing-0.5_{sm_avg}_boundary-{_gr_boundary}"), \
                f"Gelman Rubin Cut function failed with subsampling! {_gr_boundary}-{sm_avg}"
