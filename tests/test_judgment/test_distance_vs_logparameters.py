from tetres.judgment.distance_vs_logparameters import *


def test_plot_distance_posterior_comparison(thirty_taxa_chain):
    plot_distance_posteroir_comparison(thirty_taxa_chain)
    assert True


def test_plot_sos_log_comparison(thirty_taxa_chain):
    plot_sos_log_comparison(thirty_taxa_chain)
    assert True


def test_loglike_along_path(thirty_taxa_chain):
    plot_loglik_along_path(thirty_taxa_chain)
    assert True


def test_plot_log_neighbours(thirty_taxa_chain):
    for dist_type in ["kc"]:  # "rnni", "rf"]:
        plot_log_neighbours(thirty_taxa_chain, dist_type=dist_type, kc_param=0.5)
    assert True
