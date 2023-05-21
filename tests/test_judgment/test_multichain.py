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
    if platform == "linux" or platform == "linux2":
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


# def test_cen_for_each_chain(ten_taxa_multichain):
#     ten_taxa_multichain.cen_for_each_chain()
#     for chain in ten_taxa_multichain:
#         assert os.path.exists(f"{chain.working_dir}/{chain.name}_cen.tree"), "Failed to write centroid file!"


def test_clade_set_comparison(ten_taxa_multichain):
    burnin = 0.1
    try:
        os.remove(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_discrete_cc_0_1_burn-{burnin}.png")
    except FileNotFoundError:
        pass
    ten_taxa_multichain.discrete_cladesetcomparator(i=0, j=1, burnin=burnin)
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_discrete_cc_0_1_burn-{burnin}.png"), "Discrete clade set comparator failed!"


def test_clade_set_comparison_noplot(ten_taxa_multichain):
    v = ten_taxa_multichain.discrete_cladesetcomparator(i=0, j=1, plot=False)
    assert v == 99.93756243756243, "Discrete clade set comparison without plot failed!"


def test_grd_all_chains_density(ten_taxa_multichain):
    samples = 100
    try:
        os.remove(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_grd_density_full_chain_{'all' if samples == 'all' else f'subsampling_{samples}'}.pdf")
    except FileNotFoundError:
        pass
    ten_taxa_multichain.gelman_rubin_all_chains_density_plot(samples=samples)
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_grd_density_full_chain_{'all' if samples == 'all' else f'subsampling_{samples}'}.pdf"), \
        "GRD density plot failed"


def test_grd_param_choice_plot(ten_taxa_multichain):
    _gr_boundary = 0.02
    smoothing_average = "mean"
    _subsampling = False
    i = 0
    j = 1
    try:
        os.remove(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_{i}-{j}_grd_parameter_choices{f'_subsampling-{_subsampling}' if _subsampling else ''}_ess-list_{smoothing_average}_{_gr_boundary}.pdf")
    except FileNotFoundError:
        pass
    ten_taxa_multichain.gelman_rubin_parameter_choice_plot(i, j, _gr_boundary=_gr_boundary, smoothing_average=smoothing_average)
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/plots/{ten_taxa_multichain.name}_{i}-{j}_grd_parameter_choices{f'_subsampling-{_subsampling}' if _subsampling else ''}_ess-list_{smoothing_average}_{_gr_boundary}.pdf"), \
        "Failed to create the GRD parameter choice plot!"


def test_multichain_extract_cutoff(ten_taxa_multichain):
    try:
        os.remove(f"{ten_taxa_multichain.working_dir}/cutoff_files/10TaxaFix_0_10TaxaFix_1_ess-200_smoothing-0.6_mean_boundary-0.02.trees")
    except FileNotFoundError:
        pass
    try:
        os.remove(f"{ten_taxa_multichain.working_dir}/cutoff_files/10TaxaFix_0_10TaxaFix_1_ess-200_smoothing-0.6_mean_boundary-0.02.log")
    except FileNotFoundError:
        pass
    ten_taxa_multichain.extract_cutoff(i=0, j=1)
    assert os.path.exists(f"{ten_taxa_multichain.working_dir}/cutoff_files/10TaxaFix_0_10TaxaFix_1_ess-200_smoothing-0.6_mean_boundary-0.02.trees"), \
        "Failed to create cutoff file"


def test_multichain_extract_cutoff_error(ten_taxa_multichain):
    with pytest.raises(ValueError):
        # For this there was not cutoff found, ends in error
        ten_taxa_multichain.extract_cutoff(i=1, j=0)

def test_multichain_detect_burnin(ten_taxa_multichain):
        assert 10 == ten_taxa_multichain.detect_burnin(0, 1), "Burnin detection failed!"


def test_multichain_gelman_rubin_cut(ten_taxa_multichain):
    cut = ten_taxa_multichain.gelman_rubin_cut(0, 1)
    assert cut == (128, 330), "Failed GR cut!"


def test_multichain_gelman_rubin_cut_burnin(ten_taxa_multichain):
    cut = ten_taxa_multichain.gelman_rubin_cut(0, 1, burnin=50)
    assert cut == (106, 308), "Failed GR cut!"


def test_multichain_gelman_rubin_cut_burnin_subsample(ten_taxa_multichain):
    cut = ten_taxa_multichain.gelman_rubin_cut(0, 1, burnin=50, _subsampling=1)
    assert cut == (300, 704), "Failed GR cut!"
