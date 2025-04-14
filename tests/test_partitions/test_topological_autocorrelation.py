from tetres.partitions.topological_autocorrelation import topo_autocorr_plot, mean_ccd_rf


def test_topo_autocorr_plot(ten_taxa_multichain):
    topo_autocorr_plot(ten_taxa_multichain[0], ns=1000, rf=False)
    assert True


def test_mean_ccd_rf(ten_taxa_multichain):
    mean_ccd_rf(ten_taxa_multichain[0].trees)
    assert True

def test_mean_ccd_rf20(twenty_taxa_tts):
    mean_ccd_rf(twenty_taxa_tts)
    assert True
