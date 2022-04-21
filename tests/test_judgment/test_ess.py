from treeoclock.judgment.ess import coda_ess, tracerer_ess, arviz_ess


def test_coda_ess(thirty_taxa_log_data):
    assert coda_ess(thirty_taxa_log_data["posterior"],
                    chain_length=int(list(thirty_taxa_log_data["Sample"])[-1]),
                    sampling_interval=int(list(thirty_taxa_log_data["Sample"])[1])) == 1447.80470846111, "Coda ESS failed"

def test_tracerer_ess(thirty_taxa_log_data):
    assert tracerer_ess(thirty_taxa_log_data["posterior"],
                        chain_length=int(list(thirty_taxa_log_data["Sample"])[-1]),
                        sampling_interval=int(list(thirty_taxa_log_data["Sample"])[1])) == 1476.6549659664809, "Tracerer ESS failed"

def test_arviz_ess(thirty_taxa_log_data):
    assert arviz_ess(thirty_taxa_log_data["posterior"],
                     chain_length=int(list(thirty_taxa_log_data["Sample"])[-1]),
                     sampling_interval=int(list(thirty_taxa_log_data["Sample"])[1])) == 1446.1461134383974, "Arviz ESS failed"


from treeoclock.judgment.ess import _multi_dim_ess, _ess_tracerer_rsample


def test_ess_tracerer_rsample():
    _ess_tracerer_rsample()
    assert True


def test_multi_dim_ess():
    mde = _multi_dim_ess()
    assert mde == 0, "Multi dim ess failed!"
