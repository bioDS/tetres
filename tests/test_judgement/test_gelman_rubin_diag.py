from tetres.judgement._gelman_rubin_diag import gelman_rubin_cut
import os


def test_cMChain_gelman_rubin_trace(ten_taxa_multichain):
    # Testing the parameter choice plot function
    ten_taxa_multichain.gelman_rubin_parameter_choice_plot(i=0, j=1)
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_{0}-{1}_grd_parameter_choices_ess-list_mean_0.02.pdf"), \
        "Gelman Rubin Trace Plot failed!"


def test_cMChain_gelman_rubin_trace_subsampling(ten_taxa_multichain):
    # Testing the parameter choice plot function
    ten_taxa_multichain.gelman_rubin_parameter_choice_plot(i=0, j=1, _subsampling=1)
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_{0}-{1}_grd_parameter_choices_subsampling-{1}_ess-list_mean.png"), \
        "Gelman Rubin Trace Plot failed!"


def test_cMChain_gelman_rubin_plot(ten_taxa_multichain):
    # todo this test should delete the plot before recomputing it
    # todo test different sample parameters
    ten_taxa_multichain.gelman_rubin_all_chains_density_plot(samples='all')
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_grd_density_full_chain_all.pdf"), "Gelman Rubin Plot failed!"


def test_cMChain_gr_cut(ten_taxa_multichain):
    start, end = gelman_rubin_cut(ten_taxa_multichain, 0, 1, smoothing=1, ess_threshold=200, pseudo_ess_range=100, smoothing_average="mean")
    assert (start, end) == (154, 356), "Failed GR-cut"
    # Saving old cut, how does it change with the new algorithm?
    # assert (start, end) == (524, 726), "Failed GR-cut"


def test_cMChain_gr_cut_burnin(ten_taxa_multichain):
    start, end = gelman_rubin_cut(ten_taxa_multichain, 0, 1, smoothing=0.9, ess_threshold=200, pseudo_ess_range=100, smoothing_average="mean", burnin=50)
    assert (start, end) == (73, 275), "Failed GR-cut"


def test_cMChain_gr_cut_burnin_subsample(ten_taxa_multichain):
    start, end = gelman_rubin_cut(ten_taxa_multichain, 0, 1, smoothing=0.9, ess_threshold=200, pseudo_ess_range=100, smoothing_average="mean", burnin=50, _subsampling=1)
    assert (start, end) == (67, 269), "Failed GR-cut"


def test_cMChain_gr_cut_subsample(ten_taxa_multichain):
    for tolerance in [0.02, 0.01, 0.05]:
        for sm_avg in ['mean', 'median']:
            ten_taxa_multichain.gelman_rubin_cut(i=0, j=1, _overwrite=True, _subsampling=2, tolerance=tolerance, smoothing_average=sm_avg)
            assert os.path.exists(f"{ten_taxa_multichain.working_dir}/data/{ten_taxa_multichain.name}_{0}_{1}_gelman_rubin_cutoff_subsampling-2_smoothing-0.5_{sm_avg}_boundary-{tolerance}"), \
                f"Gelman Rubin Cut function failed with subsampling! {tolerance}-{sm_avg}"
