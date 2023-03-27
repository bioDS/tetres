from tetres.judgment.bic_test import bic_test


def test_bic_test(ten_taxa_multichain):
    i = 0
    j = 1
    result = bic_test(ten_taxa_multichain, i, j)
    assert result == (787.5757513291406, 790.2957635001446), "Bic test failed."
