from tetres.judgment.ess import pseudo_ess, autocorr_ess
import random


def test_pseudo_ess(thirty_taxa_MChain):
    state = random.getstate()  # get the random seed state
    random.seed(10)  # Fixing the seed to get the same result
    p_ess = int(pseudo_ess(tree_set=thirty_taxa_MChain.trees))
    random.setstate(state)  # reset the random seed to previous state
    assert p_ess == 882, "Pseudo ESS failed!"


def test_autocorr_ess(thirty_taxa_log_data):
    assert int(autocorr_ess(data_list=thirty_taxa_log_data["posterior"])) == 1478, "Autocorr ESS failed"
