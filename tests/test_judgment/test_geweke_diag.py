import pytest

from treeoclock.judgment._geweke_diag import _geweke_kind


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


# Tests for geweke_diag_summary
def test_MChain_compute_geweke_diagnostic_focal_tree_fm(thirty_taxa_MChain):
    ret = thirty_taxa_MChain.compute_geweke_focal_tree(focal_tree="FM")
    assert len(ret) == len(thirty_taxa_MChain.trees), "Geweke Focal FM failed!"


def test_MChain_compute_geweke_diagnostic_focal_tree_fm(thirty_taxa_MChain):
    ret = thirty_taxa_MChain.compute_geweke_focal_tree(focal_tree="Random")
    assert len(ret) == len(thirty_taxa_MChain.trees), "Geweke Focal Random failed!"


def test_MChain_compute_geweke_diagnostic_focal_wrongfocal(thirty_taxa_MChain):
    for ft in [12, "Bla"]:
        with pytest.raises(ValueError):
            thirty_taxa_MChain.compute_geweke_focal_tree(focal_tree=ft)


def test_MChain_compute_geweke_diagnostic_focal_wrongkind(thirty_taxa_MChain):
    for kind in [12, "Bla"]:
        with pytest.raises(ValueError):
            thirty_taxa_MChain.compute_geweke_focal_tree(kind=kind)


def test_MChain_compute_geweke_diagnostic_focal_random(thirty_taxa_MChain):
    ret = thirty_taxa_MChain.compute_geweke_focal_tree(focal_tree="Random")
    assert len(ret) == len(thirty_taxa_MChain.trees), "Geweke Focal Random failed!"


def test_MChain_compute_geweke_diagnostic_focal_FM(thirty_taxa_MChain):
    ret = thirty_taxa_MChain.compute_geweke_focal_tree(focal_tree="FM")
    assert len(ret) == len(thirty_taxa_MChain.trees), "Geweke Focal FM failed!"


# def test_MChain_compute_geweke_diagnostic_focal_centroid(thirty_taxa_MChain):
#     ret = thirty_taxa_MChain.compute_geweke_focal_tree(focal_tree="Centroid")
#     assert len(ret) == len(thirty_taxa_MChain.trees), "Geweke Focal Centroid failed!"


def test_MChain_compute_geweke_diagnostic_focal_kinds(thirty_taxa_MChain):
    for kind in _geweke_kind:
        ret = thirty_taxa_MChain.compute_geweke_focal_tree(focal_tree="Random", kind=kind)
        assert len(ret) == len(thirty_taxa_MChain.trees), f"Geweke Focal kind={kind} failed!"
