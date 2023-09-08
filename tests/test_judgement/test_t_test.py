from tetres.judgement.t_test import t_test
from tetres.summary.centroid import Centroid


def test_t_test(ten_taxa_multichain):
    cen1, sos1 = Centroid().compute_centroid(ten_taxa_multichain[0].trees)
    cen2, sos2 = Centroid().compute_centroid(ten_taxa_multichain[1].trees)
    test_stats = t_test(trees1=ten_taxa_multichain[0].trees, trees2=ten_taxa_multichain[1].trees, centroid1=cen1, centroid2=cen2)
    assert test_stats == (3.2981027725073995, 0.0009904691869427231), "Something went wrong."
