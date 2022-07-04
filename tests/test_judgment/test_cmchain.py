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


def test_cMChain_pwd_matrix_rf(ten_taxa_cMChain):
    errors = []
    dm0 = ten_taxa_cMChain.pwd_matrix(0, rf=True)
    if dm0.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    dm1 = ten_taxa_cMChain.pwd_matrix(1, rf=True)
    if dm1.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    dm01 = ten_taxa_cMChain.pwd_matrix(0, 1, rf=True)
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


def test_cMChain_gress(ten_taxa_cMChain):
    ten_taxa_cMChain.gelman_rubin_trace_ess_plot(i=0, j=1)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_{0}-{1}_gress_arviz.png"), \
        "Gelman Rubin Trace Plot with ess failed!"


def test_cMChain_cladesetcomparator(ten_taxa_cMChain):
    ten_taxa_cMChain.cladesetcomparator(beast_applauncher="/home/lars/.local/share/beast/bin/applauncher")
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/cc.png"), "Cladeset comparator failed!"


def test_cMChain_ess_stripplot(ten_taxa_cMChain):
    for method in ["tracerer", "coda", "arviz"]:
        ten_taxa_cMChain.ess_stripplot(ess_method=method)
        assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/ess_{method}_comparison.png"), "ESS stripplot failed!"


# todo test temporary
def test_cMChain_spectral_cluster_all(ten_taxa_cMChain):
    ten_taxa_cMChain.spectral_cluster_all(n_clus=3)
    assert True


def test_cMChain_split_all_trees(ten_taxa_cMChain):
    ten_taxa_cMChain.split_all_trees(n_clus=3)
    # todo this test should first remove the files and then run the test
    assert True


def test_pwd_matrix_all(ten_taxa_cMChain):
    matrix = ten_taxa_cMChain.pwd_matrix_all()
    assert matrix.shape == (3003, 3003), "PWD matrix all failed!"


def test_pwd_matrix_all_rf(ten_taxa_cMChain):
    matrix = ten_taxa_cMChain.pwd_matrix_all(rf=True)
    assert matrix.shape == (3003, 3003), "PWD matrix RF all failed!"
