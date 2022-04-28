from treeoclock.judgment._plots import _ess_trace_plot

import pandas as pd
import numpy as np


def test_ess_trace_plot_pseudo_ess(thirty_taxa_MChain):
    data = []
    for i in range(5, len(thirty_taxa_MChain.trees)):
        data.append(["Pseudo_ess", thirty_taxa_MChain.get_pseudo_ess(upper_i=i), i])
    data = pd.DataFrame(data, columns=["Ess_key", "Ess_value", "Upper_i"])
    assert _ess_trace_plot(data) == 0, "_ess_trace_plot failed!"


def test_ess_trace_plot_logfile_ess(thirty_taxa_MChain):
    for ess_method in ['tracerer', 'arviz', 'coda']:
        data = []
        ess_key = thirty_taxa_MChain.get_key_names()
        for i in range(5, thirty_taxa_MChain.log_data.shape[0]):  # Starting at sample 5 as it is not useful to look at less samples
            data.extend([[key, thirty_taxa_MChain.get_ess(ess_key=key, ess_method=ess_method, upper_i=i), i] for key in ess_key])
        data = pd.DataFrame(data, columns=["Ess_key", "Ess_value", "Upper_i"])
        assert _ess_trace_plot(data) == 0, f"_ess_trace_plot using {ess_method} failed!"


def test_ess_trace_plot_with_rnni_variance(thirty_taxa_MChain):
    for ess_method in ['tracerer', 'arviz', 'coda']:
        data = []
        thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="tree")
        thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="FM")
        ess_key = thirty_taxa_MChain.get_key_names()
        for i in range(5, thirty_taxa_MChain.log_data.shape[0]):  # Starting at sample 5 as it is not useful to look at less samples
                data.extend(
                    [[key, thirty_taxa_MChain.get_ess(ess_key=key, ess_method=ess_method, upper_i=i), i] for key in ess_key
                     if not np.isnan(thirty_taxa_MChain.log_data[key][i])])
        data = pd.DataFrame(data, columns=["Ess_key", "Ess_value", "Upper_i"])
        assert _ess_trace_plot(data) == 0, f"RNNI variance Plot went wrong!"


def test_ess_trace_plot_with_rnni_distance(thirty_taxa_MChain):
    for ess_method in ['tracerer', 'arviz', 'coda']:
        data = []
        for average in ["mean", "median", "median_ad", "mean_ad"]:
            thirty_taxa_MChain.compute_new_tree_distance_log(average=average)
        ess_key = thirty_taxa_MChain.get_key_names()
        for i in range(5, thirty_taxa_MChain.log_data.shape[0]):  # Starting at sample 5 as it is not useful to look at less samples
                data.extend(
                    [[key, thirty_taxa_MChain.get_ess(ess_key=key, ess_method=ess_method, upper_i=i), i] for key in ess_key
                     if not np.isnan(thirty_taxa_MChain.log_data[key][i])])
        data = pd.DataFrame(data, columns=["Ess_key", "Ess_value", "Upper_i"])
        assert _ess_trace_plot(data) == 0, f"RNNI new tree distance log plot went wrong!"


def test_ess_trace_plot_with_rnni_summary_distance(thirty_taxa_MChain):
    for ess_method in ['tracerer', 'arviz', 'coda']:
        data = []
        thirty_taxa_MChain.compute_new_tree_summary_distance_log(summary="FM")
        ess_key = thirty_taxa_MChain.get_key_names()
        for i in range(5, thirty_taxa_MChain.log_data.shape[0]):  # Starting at sample 5 as it is not useful to look at less samples
                data.extend(
                    [[key, thirty_taxa_MChain.get_ess(ess_key=key, ess_method=ess_method, upper_i=i), i] for key in ess_key
                     if not np.isnan(thirty_taxa_MChain.log_data[key][i])])
        data = pd.DataFrame(data, columns=["Ess_key", "Ess_value", "Upper_i"])
        assert _ess_trace_plot(data) == 0, f"RNNI distance summary new tree went wrong!"
