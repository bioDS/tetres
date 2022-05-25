import pytest
from pathlib import Path
from treeoclock.judgment.mchain import coupled_MChains


def test_coupledMChain():
    assert coupled_MChains(m_MChains=1,
                           trees=["30Taxa.trees"],
                           log_files=["30Taxa_beast2.log"],
                           working_dir= f"{Path(__file__).parent.parent.absolute()}/data"
                           ), "Construction coupledMChain 1 with files failed!"


def test_coupledMChain():
    assert coupled_MChains(m_MChains=2,
                           trees=["chain0.trees", "chain0_1.trees"],
                           log_files=["chain0.log", "chain0_1.log"],
                           working_dir= f"{Path(__file__).parent.parent.absolute()}/data/cMChain"
                           ), "Construction coupledMChain 2 with files failed!"


def test_cMChain_gelman_rubin_plot_single_chain():
    with pytest.raises(ValueError):
        coupled_chains = coupled_MChains(m_MChains=1,
                        trees=["30Taxa.trees"],
                        log_files=["30Taxa_beast2.log"],
                        working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
                        )
        coupled_chains.gelman_rubin_like_diagnostic_plot()


def test_cMChain_gelman_rubin_plot(ten_taxa_cMChain):
    assert ten_taxa_cMChain.gelman_rubin_like_diagnostic_plot(), "Gelman Rubin Plot failed!"



