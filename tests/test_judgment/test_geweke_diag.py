import pytest


# Tests for geweke_diag_distances
def test_MChain_compute_geweke_distances(thirty_taxa_MChain):
    ret = thirty_taxa_MChain.compute_geweke_distances()
    assert len(ret) == len(thirty_taxa_MChain.trees), "Geweke_distances failed!"


def test_MChain_compute_geweke_distances_diffintervals(thirty_taxa_MChain):
    ret = thirty_taxa_MChain.compute_geweke_distances(first_range=[0, 0.25], last_percent=0.5)
    assert len(ret) == len(thirty_taxa_MChain.trees), "Geweke_distances failed!"


def test_MChain_compute_geweke_distances_wrongintervals(thirty_taxa_MChain):
    for first, last in [([1, 2, 3], 0.4), ([0.2, 0.1], 0.4), ([0.1, 0.1], 0.4), ([1, 1], 0.4), ([0.1, 0.2], 1), ([0.1, 0.35], 0.7)]:
        with pytest.raises(ValueError):
            thirty_taxa_MChain.compute_geweke_distances(first_range=first, last_percent=last)

# Tests for geweke_diag_



