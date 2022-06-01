import pytest
from pathlib import Path
from treeoclock.judgment.mchain import coupled_MChains
import os


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


def test_cMChain_pwd_matrix(ten_taxa_cMChain):
    errors = []
    dm0 = ten_taxa_cMChain.pwd_matrix(0)
    if dm0.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    dm1 = ten_taxa_cMChain.pwd_matrix(1)
    if dm1.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    dm01 = ten_taxa_cMChain.pwd_matrix(0, 1)
    if dm01.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    assert errors == [], f"PWD matrix went wrong! {errors}"


def test_cMChain_pwdm_errors(ten_taxa_cMChain):
    with pytest.raises(ValueError):
        for i in ["HI", 1.0, 2]:
            ten_taxa_cMChain.pwd_matrix(i)
            ten_taxa_cMChain.pwd_matrix(0, i)
        ten_taxa_cMChain.pw_matrix(0, 0)


def test_cMChain_gelman_rubin_plot(ten_taxa_cMChain):
    ten_taxa_cMChain.gelman_rubin_like_diagnostic_plot()
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_grd_plot.png"), "Gelman Rubin Plot failed!"


def test_cMChain_gelman_rubin_trace(ten_taxa_cMChain):
    ten_taxa_cMChain.gelman_rubin_trace_plot(i=0, j=1)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_grd_singlevalue_evaluation.png"), \
        "Gelman Rubin Trace Plot failed!"
