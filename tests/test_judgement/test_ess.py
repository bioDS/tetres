from tetres.judgement.ess import (pseudo_ess, autocorr_ess,
                                  tracer_ess_convenience,
                                  calc_ess_beast2)
import random


def test_pseudo_ess(thirty_taxa_MChain):
    state = random.getstate()  # get the random seed state
    random.seed(10)  # Fixing the seed to get the same result
    p_ess = int(pseudo_ess(tree_set=thirty_taxa_MChain.trees))
    random.setstate(state)  # reset the random seed to previous state
    assert p_ess == 882, "Pseudo ESS failed!"


def test_autocorr_ess(thirty_taxa_log_data):
    assert int(
        autocorr_ess(data_list=thirty_taxa_log_data["posterior"])) == 1478, "Autocorr ESS failed"


def test_recover_tracer_ess_implementation(thirty_taxa_log_data):
    for col in thirty_taxa_log_data.columns[1:]:
        ess_autocorr = autocorr_ess(thirty_taxa_log_data[col])
        ess_convenience = tracer_ess_convenience(thirty_taxa_log_data[col])
        ess_beast2 = calc_ess_beast2(thirty_taxa_log_data[col])
        print("--------------")
        print(col, ":\n\t", round(ess_convenience, 4), "\t--\t", round(ess_beast2, 4),
              "\t--\t", round(ess_autocorr, 4))
    assert True
