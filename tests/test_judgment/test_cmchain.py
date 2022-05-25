import pytest
from pathlib import Path
from treeoclock.judgment.mchain import coupled_MChains


def test_coupledMChain():
    assert coupled_MChains(m_MChains=1,
                           trees=["30Taxa.trees"],
                           log_files=["30Taxa_beast2.log"],
                           working_dir= f"{Path(__file__).parent.parent.absolute()}/data"
                           ), "Construction coupledMChain with files failed!"


def test_cMChain_gelman_rubin_plot_single_chain():
    with pytest.raises(ValueError):
        coupled_chains = coupled_MChains(m_MChains=1,
                        trees=["30Taxa.trees"],
                        log_files=["30Taxa_beast2.log"],
                        working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
                        )
        coupled_chains.gelman_rubin_like_diagnostic_plot()


