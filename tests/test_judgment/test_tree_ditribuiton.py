import pytest
import os
from treeoclock.judgment.tree_distribution import get_coefficients, _get_approx_orbits, plot_tree_density_distribution

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
    plot_tree_density_distribution(Mchain=ten_taxa_cMChain, given_x=10)
    assert True
