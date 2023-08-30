from tetres.partitions.topological_autocorrelation import *


def test_topo_autocorr_plot(ten_taxa_multichain):
    topo_autocorr_plot(ten_taxa_multichain[0], ns=1000, rf=False)
    assert True
