from tetres.partitions.HD_traces import *


def test_hd_trace(ten_taxa_multichain):
    data = entropy_trace_data(ten_taxa_multichain)
    assert True

def test_hd_ranked_trace(ten_taxa_multichain):
    data = entropy_trace_ranked_data(ten_taxa_multichain)
    assert True
