from treeoclock.judgment.ess import coda_ess, tracerer_ess, arviz_ess, pseudo_ess
import random


def test_coda_ess(thirty_taxa_log_data):
    assert coda_ess(data_list=thirty_taxa_log_data["posterior"],
                    chain_length=int(list(thirty_taxa_log_data["Sample"])[-1]),
                    sampling_interval=int(
                        list(thirty_taxa_log_data["Sample"])[1])) == 1447.80470846111, "Coda ESS failed"


def test_tracerer_ess(thirty_taxa_log_data):
    assert tracerer_ess(data_list=thirty_taxa_log_data["posterior"],
                        chain_length=int(list(thirty_taxa_log_data["Sample"])[-1]),
                        sampling_interval=int(
                            list(thirty_taxa_log_data["Sample"])[1])) == 1476.6549659664809, "Tracerer ESS failed"


def test_arviz_ess(thirty_taxa_log_data):
    assert arviz_ess(data_list=thirty_taxa_log_data["posterior"],
                     chain_length=int(list(thirty_taxa_log_data["Sample"])[-1]),
                     sampling_interval=int(
                         list(thirty_taxa_log_data["Sample"])[1])) == 1446.1461134383974, "Arviz ESS failed"


def test_pseudo_ess(thirty_taxa_MChain):
    state = random.getstate()  # get the random seed state
    random.seed(10)  # Fixing the seed to get the same result
    p_ess = []
    for method in ["tracerer", "arviz", "coda"]:
        p_ess.append(int(pseudo_ess(ess_method="tracerer", tree_set=thirty_taxa_MChain.trees, chain_length=thirty_taxa_MChain.chain_length,
                           sampling_interval=thirty_taxa_MChain.chain_length / (len(thirty_taxa_MChain.trees) - 1))))
    random.setstate(state)  # reset the random seed to previous state
    assert p_ess == [865, 855, 866], "Pseudo ESS failed!"


from treeoclock.judgment.ess import _multi_dim_ess, _ess_tracerer_rsample


def test_ess_tracerer_rsample():
    _ess_tracerer_rsample()
    assert True


def test_multi_dim_ess():
    mde = _multi_dim_ess()
    assert mde == 0, "Multi dim ess failed!"
