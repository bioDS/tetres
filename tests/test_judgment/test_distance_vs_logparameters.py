from treeoclock.judgment.distance_vs_logparameters import plot_distance_posteroir_comparison, plot_sos_log_comparison

def test_plot_distance_posterior_comparison(thirty_taxa_MChain):
    plot_distance_posteroir_comparison(thirty_taxa_MChain)
    assert True

def test_plot_sos_log_comparison(thirty_taxa_MChain):
    plot_sos_log_comparison(thirty_taxa_MChain)
    assert True
