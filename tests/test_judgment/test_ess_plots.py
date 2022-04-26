from treeoclock.judgment._ess_plots import _ess_trace_plot

import pandas as pd

# todo this is just temporary
def test_ess_trace_plot_pseudo_ess(thirty_taxa_MChain):
    data = []
    for i in range(5, len(thirty_taxa_MChain.trees)):
        data.append(["Pseudo_ess", thirty_taxa_MChain.get_pseudo_ess(upper_i=i), i])
    data = pd.DataFrame(data, columns=["Ess_key", "Ess_value", "Upper_i"])
    assert _ess_trace_plot(data) == 0, "_ess_trace_plot failed!"

