from tetres.judgment.bic_test import bic_test


def test_bic_test(ten_taxa_multichain):
    result = bic_test(ten_taxa_multichain, i, j)
    assert result == 0
