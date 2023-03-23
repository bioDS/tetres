import pytest
import os
from tetres.judgment.tree_distribution import get_coefficients,\
    _get_approx_orbits, plot_tree_density_distribution, \
    plot_CCD_vs_centroid_distance, get_sample_treespace_coverage, \
    get_corr_p_coverage

def test_sympy_product():
    o = get_coefficients(5)
    o2 = get_coefficients(6)
    assert (sum(o), sum(o2)) == (180, 2700), "Orbit polynom for caterpillar failed!"


def test_get_approx_orbits():
    # for 5 taxa there are 180 trees in total
    o = _get_approx_orbits(5)
    o2 = _get_approx_orbits(6)
    assert (sum(o), sum(o2)) == (180, 2700), "Orbit polynom for caterpillar failed!"


def test_plot_tree_density_distribution(ten_taxa_cMChain):
    plot_tree_density_distribution(Mchain=ten_taxa_cMChain)
    assert True


def test_get_sample_treespace_coverage(ten_taxa_cMChain):
    c = get_sample_treespace_coverage(Mchain=ten_taxa_cMChain)
    assert c == 0.00035984908828859896, "Treespace covereage failed!"


def test_plot_CCD_vs_centroid_distance(ten_taxa_cMChain):
    plot_CCD_vs_centroid_distance(ten_taxa_cMChain)
    assert True


def test_get_corr_p_coverage(ten_taxa_cMChain):
    pr, pr_f, pr_t, coverage = get_corr_p_coverage(Mchain=ten_taxa_cMChain)
    assert True
