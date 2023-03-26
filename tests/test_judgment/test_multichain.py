import pytest
from pathlib import Path
from tetres.judgment.multichain import MultiChain
import os
from sys import platform


def test_multichain():
    assert MultiChain(m_chains=1,
                           trees=["30Taxa.trees"],
                           log_files=["30Taxa_beast2.log"],
                           working_dir= f"{Path(__file__).parent.parent.absolute()}/data"
                           ), "Construction coupledMChain 1 with files failed!"


def test_multichain():
    assert MultiChain(m_chains=2,
                           trees=["chain0.trees", "chain0_1.trees"],
                           log_files=["chain0.log", "chain0_1.log"],
                           working_dir= f"{Path(__file__).parent.parent.absolute()}/data/cMChain"
                           ), "Construction coupledMChain 2 with files failed!"


def test_diff_map():
    with pytest.raises(ValueError):
        MultiChain(m_chains=2,
                trees=["chain0.trees", "chain0_1_diff_map.trees"],
                log_files=["chain0.log", "chain0_1.log"],
                working_dir= f"{Path(__file__).parent.parent.absolute()}/data/cMChain", name="Diff"
                )


def test_multichain_get_ess_start_end(ten_taxa_multichain):
    ess = ten_taxa_multichain[0].get_ess(ess_key="posterior", lower_i=10, upper_i=110)
    assert int(ess) == 67, "Get ESS with start end failed"


def test_multichain_pwd_matrix(ten_taxa_multichain):
    errors = []
    dm0 = ten_taxa_multichain.pwd_matrix(0)
    if dm0.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    dm1 = ten_taxa_multichain.pwd_matrix(1)
    if dm1.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    dm01 = ten_taxa_multichain.pwd_matrix(0, 1)
    if dm01.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    assert errors == [], f"PWD matrix went wrong! {errors}"


def test_multichain_pwd_matrix_rf(ten_taxa_multichain):
    errors = []
    dm0 = ten_taxa_multichain.pwd_matrix(0, rf=True)
    if dm0.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    dm1 = ten_taxa_multichain.pwd_matrix(1, rf=True)
    if dm1.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    dm01 = ten_taxa_multichain.pwd_matrix(0, 1, rf=True)
    if dm01.shape != (1001, 1001):
        errors.append("dm0 wrong shape!")
    assert errors == [], f"PWD matrix went wrong! {errors}"


def test_multichain_pwdm_errors(ten_taxa_multichain):
    with pytest.raises(ValueError):
        for i in ["HI", 1.0, 2]:
            ten_taxa_multichain.pwd_matrix(i)
            ten_taxa_multichain.pwd_matrix(0, i)
        ten_taxa_multichain.pw_matrix(0, 0)


def test_multichain_cladesetcomparator(ten_taxa_multichain):
    ten_taxa_multichain.cladesetcomparator(beast_applauncher="/Applications/BEAST 2.6.3/bin/applauncher")
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/plots/10TaxaFix_clade_comp.png"), "Cladeset comparator failed!"


def test_pwd_matrix_all(ten_taxa_multichain):
    if platform == "linux" or platform == "linux2":
        for rf in [True, False]:
            matrix = ten_taxa_multichain.pwd_matrix_all(rf=rf)
            assert matrix.shape == (3003, 3003), "PWD matrix all failed!"
    else:
        matrix = ten_taxa_multichain.pwd_matrix_all(rf=False)
        assert matrix.shape == (3003, 3003), "PWD matrix all failed!"


def test_pwd_matrix_all_mac_warning(ten_taxa_multichain):
    if platform == "darwin" or platform == "win32":
        with pytest.raises(ValueError):
            ten_taxa_multichain.pwd_matrix_all(rf=True)


def test_cen_for_each_chain(ten_taxa_multichain):
    ten_taxa_multichain.cen_for_each_chain()
    for chain in ten_taxa_multichain:
        assert os.path.exists(f"{chain.working_dir}/{chain.name}_cen.tree"), "Failed to write centroid file!"


def test_clade_set_comparison(ten_taxa_multichain):
    burnin = 0.1
    try:
        os.remove(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_discrete_cc_0_1_burn-{burnin}.png")
    except FileNotFoundError:
        pass
    ten_taxa_multichain.clade_set_comparison(i=0, j=1, burnin=burnin)
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_discrete_cc_0_1_burn-{burnin}.png"), "Discrete clade set comparator failed!"


def test_clade_set_comparison_noplot(ten_taxa_multichain):
    v = ten_taxa_multichain.clade_set_comparison(i=0, j=1, plot=False)
    assert v == 99.93756243756243, "Discrete clade set comparison without plot failed!"


def test_get_best_bic_cluster(ten_taxa_multichain):
    best_bic = ten_taxa_multichain[0].get_best_bic_cluster()
    assert best_bic == 5, "Get best bic Mchain function failed"


def test_get_best_silhouette_cluster(ten_taxa_multichain):
    best_sil = ten_taxa_multichain[0].get_best_silhouette_cluster()
    assert best_sil == 5, "Get best silhouette cluster evaluation failed"
