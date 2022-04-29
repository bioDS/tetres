import pytest
import random
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from treeoclock.judgment.mchain import MChain
import treeoclock.judgment.ess as ess


def test_MChain_construction_files():
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction files failed"


def test_MChain_construction_TimeTreeSets(thirty_taxa_tts, thirty_taxa_centroid_tts):
    assert MChain(
        trees=thirty_taxa_tts,
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=thirty_taxa_centroid_tts,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction TTS failed"


def test_MChain_construction_mixed1(thirty_taxa_tts):
    assert MChain(
        trees=thirty_taxa_tts,
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction TTS trees, centroid file failed"


def test_MChain_construciton_mixed2(thirty_taxa_centroid_tts):
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=thirty_taxa_centroid_tts,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction trees file, TTS centroid failed"


def test_MChain_construciton_summary_None():
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=None,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction with summary=None failed"


def test_MChain_construciton_workingdir_None():
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
        working_dir=None
    ), "Construction with working_dir=None failed"


def test_MChain_construciton_workingdir_summary_None():
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=None,
        working_dir=None
    ), "Construction with working_dir=None failed"


# testing construction exceptions
def test_MChain_wrong_trees():
    with pytest.raises(ValueError):
        MChain(
            trees=20,
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_wrong_logfile_path():
    with pytest.raises(FileNotFoundError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_wrong_logfile_type():
    with pytest.raises(ValueError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=20,
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_wrong_summary():
    with pytest.raises(ValueError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=20,
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_wrong_workingdir_path():
    with pytest.raises(NotADirectoryError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/work/"
        )


def test_MChain_wrong_workingdir_type():
    with pytest.raises(ValueError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=20
        )


def test_MChain_wrong_treemaps():
    with pytest.raises(ValueError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/12Taxa.trees",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_get_key_names(thirty_taxa_MChain):
    assert thirty_taxa_MChain.get_key_names() == ['posterior',
                                                  'likelihood',
                                                  'prior',
                                                  'treeLikelihood',
                                                  'TreeHeight',
                                                  'popSize',
                                                  'CoalescentConstant'], \
        "Get_key_name functions failed!"


def test_MChain_get_ess_tracerer(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="tracerer"))
    assert result == [1476.6549659664809,
                      1203.6766888338857,
                      1668.228049029279,
                      1203.6766888338857,
                      1758.881844170911,
                      1650.888969360209,
                      1677.680147437517], \
        "ESS tracerer value calculation with MChain failed!"


def test_MChain_get_ess_coda(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="coda"))
    assert result == [1447.80470846111,
                      1398.5896384868036,
                      1697.2455114765976,
                      1398.5896384868036,
                      1819.6601044141469,
                      1635.0719443144444,
                      1707.950568955924], \
        "ESS coda value calculation with MChain failed!"


def test_MChain_get_ess_arviz(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="arviz"))
    assert result == [1446.1461134383974,
                      1388.2577408057957,
                      1680.2897556479547,
                      1388.2577408057957,
                      1830.3012635184548,
                      1661.0004487633273,
                      1683.7924322905303], \
        "ESS arviz value calculation with MChain failed!"


def test_MChain_get_ess_partial_tracerer(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="tracerer", upper_i=200, lower_i=100))
    assert result == [100.0, 90.79172611567039, 100.0, 90.79172611567039, 100.0, 100.0, 100.0], \
        "ESS partial tracerer value calculation with MChain failed!"


def test_MChain_get_ess_partial_coda(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="coda", upper_i=200, lower_i=100))
    assert result == [99.99999999999999,
                      100.00000000000001,
                      100.00000000000001,
                      100.00000000000001,
                      99.99999999999996,
                      132.39259058042072,
                      100.0], \
        "ESS partial coda value calculation with MChain failed!"


def test_MChain_get_ess_partial_arviz(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="arviz", upper_i=200, lower_i=100))
    assert result == [122.731291995835,
                      95.5331836361955,
                      131.99641052649312,
                      95.5331836361955,
                      93.74009808720129,
                      152.44302830034806,
                      130.80886336654285], \
        "ESS partial arviz value calculation with MChain failed!"


def test_MChain_get_ess_partial_wronglowertype(thirty_taxa_MChain):
    with pytest.raises(ValueError):
        thirty_taxa_MChain.get_ess(ess_key="posterior", ess_method="arviz", lower_i=10.2)


def test_MChain_get_ess_partial_wronguppertype(thirty_taxa_MChain):
    with pytest.raises(ValueError):
        thirty_taxa_MChain.get_ess(ess_key="posterior", ess_method="arviz", upper_i=10.2)


def test_MChain_get_ess_partial_uppertoobig(thirty_taxa_MChain):
    with pytest.raises(IndexError):
        thirty_taxa_MChain.get_ess(ess_key="posterior", ess_method="arviz", upper_i=3000)


def test_MChain_get_ess_partial_wrongway(thirty_taxa_MChain):
    with pytest.raises(ValueError):
        thirty_taxa_MChain.get_ess(ess_key="posterior", ess_method="arviz", lower_i=10, upper_i=5)


def test_MChain_get_ess_partial_neglower(thirty_taxa_MChain):
    with pytest.raises(ValueError):
        thirty_taxa_MChain.get_ess(ess_key="posterior", ess_method="arviz", lower_i=-10)


def test_MChain_get_ess_partial_negupper(thirty_taxa_MChain):
    with pytest.raises(ValueError):
        thirty_taxa_MChain.get_ess(ess_key="posterior", ess_method="arviz", upper_i=-10)


def test_MChain_get_ess_partial_correctness_arviz(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="arviz"))

    result_partial = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result_partial.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="arviz", lower_i=0,
                                                         upper_i=thirty_taxa_MChain.log_data.shape[0]))
    assert result == result_partial, \
        "ESS partial returns wrong ESS!"


def test_MChain_get_ess_partial_correctness_coda(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="coda"))

    result_partial = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result_partial.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="coda", lower_i=0,
                                                         upper_i=thirty_taxa_MChain.log_data.shape[0]))
    assert result == result_partial, \
        "ESS partial returns wrong ESS!"


def test_MChain_get_ess_partial_correctness_tracerer(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="tracerer"))

    result_partial = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result_partial.append(thirty_taxa_MChain.get_ess(ess_key=ess_key, ess_method="tracerer", lower_i=0,
                                                         upper_i=thirty_taxa_MChain.log_data.shape[0]))
    assert result == result_partial, \
        "ESS partial returns wrong ESS!"


def test_MChain_get_ess_trace_plot_default(thirty_taxa_MChain, monkeypatch):
    monkeypatch.setattr(plt, 'show', lambda: None)  # surpress plt.show() for plot, only temporary
    assert thirty_taxa_MChain.get_ess_trace_plot() == 0, "ESS trace plot funciton failed"


def test_MChain_get_ess_trace_plot_one_key(thirty_taxa_MChain, monkeypatch):
    monkeypatch.setattr(plt, 'show', lambda: None)  # surpress plt.show() for plot, only temporary
    assert thirty_taxa_MChain.get_ess_trace_plot(ess_key="posterior") == 0, "ESS trace plot funciton failed"


def test_MChain_get_ess_trace_plot_list_key(thirty_taxa_MChain, monkeypatch):
    monkeypatch.setattr(plt, 'show', lambda: None)  # surpress plt.show() for plot, only temporary
    assert thirty_taxa_MChain.get_ess_trace_plot(ess_key=["posterior", "likelihood"]) == 0, "ESS trace plot funciton failed"


def test_MChain_get_ess_trace_plot_wrong_list(thirty_taxa_MChain, monkeypatch):
    with pytest.raises(ValueError):
        thirty_taxa_MChain.get_ess_trace_plot(ess_key=["posterior", "Loglik"])


def test_MChain_get_ess_trace_plot_wrong_string(thirty_taxa_MChain, monkeypatch):
    with pytest.raises(ValueError):
        thirty_taxa_MChain.get_ess_trace_plot(ess_key="Posterior")


def test_MChain_get_ess_trace_plot_wrong_list2(thirty_taxa_MChain, monkeypatch):
    with pytest.raises(ValueError):
        thirty_taxa_MChain.get_ess_trace_plot(ess_key=["Posterior", 2, 2.003])


def test_MChain_get_trace_plot_new_tree_distance(thirty_taxa_MChain, monkeypatch):
    for average in ["mean", "median", "median_ad", "mean_ad"]:
        thirty_taxa_MChain.compute_new_tree_distance_log(average=average)
    monkeypatch.setattr(plt, 'show', lambda: None)  # surpress plt.show() for plot, only temporary
    for average in ["mean", "median", "median_ad", "mean_ad"]:
        assert thirty_taxa_MChain.get_trace_plot(f"Distance_{average}") == 0, \
            f"Trace Plot for Distance_{average} went wrong!"


def test_MChain_get_trace_plot_wrong_key(thirty_taxa_MChain):
    with pytest.raises(KeyError):
        thirty_taxa_MChain.get_trace_plot("Posterior")


def test_MChain_pseudo_ess(thirty_taxa_MChain):
    state = random.getstate()  # get the random seed state
    random.seed(10)  # Fixing the seed to get the same result
    p_ess = thirty_taxa_MChain.get_pseudo_ess()
    random.setstate(state)  # reset the random seed to previous state
    assert int(p_ess) == 921, "Get Pseudo ESS for MChain failed!"


def test_MChain_compute_distance_ess_log_arviz(thirty_taxa_MChain):
    ess_values = []
    ess_method = "arviz"
    for average in ["mean", "median", "median_ad", "mean_ad"]:
        thirty_taxa_MChain.compute_new_tree_distance_log(average=average)
        ess_values.append(int(thirty_taxa_MChain.get_ess(ess_key=f"Distance_{average}", ess_method=ess_method)))
    assert ess_values == [843, 853, 638, 763], "Compute distance ess log values failed!"


def test_MChain_compute_distance_ess_log_coda(thirty_taxa_MChain):
    ess_values = []
    ess_method = "coda"
    for average in ["mean", "median", "median_ad", "mean_ad"]:
        thirty_taxa_MChain.compute_new_tree_distance_log(average=average)
        ess_values.append(int(thirty_taxa_MChain.get_ess(ess_key=f"Distance_{average}", ess_method=ess_method)))
    assert ess_values == [997, 1013, 252, 311], "Compute distance ess log values failed!"


def test_MChain_compute_distance_ess_log_tracerer(thirty_taxa_MChain):
    ess_values = []
    ess_method = "tracerer"
    for average in ["mean", "median", "median_ad", "mean_ad"]:
        thirty_taxa_MChain.compute_new_tree_distance_log(average=average)
        ess_values.append(int(thirty_taxa_MChain.get_ess(ess_key=f"Distance_{average}", ess_method=ess_method)))
    assert ess_values == [972, 978, 247, 302], "Compute distance ess log values failed!"


# todo missing more tests for other summary options
def test_MChain_compute_distance_new_tree_summary_log_tracerer(thirty_taxa_MChain):
    ess_method = "tracerer"
    summary = "FM"
    thirty_taxa_MChain.compute_new_tree_summary_distance_log(summary=summary)
    ess_value = int(thirty_taxa_MChain.get_ess(ess_key=f"Distance_{summary}", ess_method=ess_method))
    assert ess_value == 911, "Compute distance ess log values failed!"


def test_MChain_compute_distance_new_tree_summary_log_arviz(thirty_taxa_MChain):
    ess_method = "arviz"
    summary = "FM"
    thirty_taxa_MChain.compute_new_tree_summary_distance_log(summary=summary)
    ess_value = int(thirty_taxa_MChain.get_ess(ess_key=f"Distance_{summary}", ess_method=ess_method))
    assert ess_value == 636, "Compute distance ess log values failed!"


def test_MChain_compute_distance_new_tree_summary_log_coda(thirty_taxa_MChain):
    ess_method = "coda"
    summary = "FM"
    thirty_taxa_MChain.compute_new_tree_summary_distance_log(summary=summary)
    ess_value = int(thirty_taxa_MChain.get_ess(ess_key=f"Distance_{summary}", ess_method=ess_method))
    assert ess_value == 551, "Compute distance ess log values failed!"

# todo
#  Missing tests for the compute_rnni_variance_log
#  Missing tests for the add_new_loglist


# todo these two tests are deprecated!

def test_MChain_compute_new_log_data_RNNIVariance(thirty_taxa_MChain):
    log_list = thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="tree")
    assert [type(log_list), len(log_list)] == [list, len(thirty_taxa_MChain.trees)], "New log list computation failed!"


def test_MChain_compute_new_log_data_RNNIVariance_ess(thirty_taxa_MChain):
    state = random.getstate()  # get the random seed state
    random.seed(10)  # Fixing the seed to get the same result
    log_list = thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="tree")
    random.setstate(state)  # reset the random seed to previous state
    ess_values = []
    for ess_method in ["arviz", "coda", "tracerer"]:
        ess_values.append(int(getattr(ess, f"{ess_method}_ess")(
            data_list=log_list, chain_length=thirty_taxa_MChain.chain_length,
            sampling_interval=thirty_taxa_MChain.chain_length/(len(log_list)-1))))
    assert ess_values == [764, 779, 611], "ESS rnniVariance list failed"


# todo this is a very lengthy test on this dataset, maybe i need a smaller test dataset for these things
def test_MChain_write_log_file(thirty_taxa_MChain):
    thirty_taxa_MChain.compute_new_tree_summary_distance_log(summary="FM")
    thirty_taxa_MChain.compute_new_tree_summary_distance_log(summary="FM", norm=True)
    for average in ["mean", "median", "median_ad", "mean_ad"]:
        thirty_taxa_MChain.compute_new_tree_distance_log(average=average)
        thirty_taxa_MChain.compute_new_tree_distance_log(average=average, norm=True)
    thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="tree")
    thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="tree", norm=True)
    thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="FM")
    thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="FM", norm=True)
    thirty_taxa_MChain.compute_ess_traces()
    thirty_taxa_MChain.write_log_file(f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_pd.log")


def test_MChain_copute_ess_traces(thirty_taxa_MChain):
    thirty_taxa_MChain.compute_new_tree_distance_log(average="mean_ad", norm=True)
    thirty_taxa_MChain.compute_rnni_variance_log(focal_tree_type="tree", norm=True)
    assert thirty_taxa_MChain.compute_ess_traces() == 0, "Computing ess traces failed!"


