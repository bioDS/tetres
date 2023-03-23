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


def test_cMChain_get_ess_start_end(ten_taxa_cMChain):
    ess = ten_taxa_cMChain[0].get_ess(ess_key="posterior", lower_i=10, upper_i=110)
    assert int(ess) == 66, "Get ESS with start end failed"


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


def test_cMChain_cladesetcomparator(ten_taxa_cMChain):
    ten_taxa_cMChain.cladesetcomparator(beast_applauncher="/home/lars/.local/share/beast/bin/applauncher")
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/cc.png"), "Cladeset comparator failed!"


def test_pwd_matrix_all(ten_taxa_cMChain):
    for rf in [True, False]:
        matrix = ten_taxa_cMChain.pwd_matrix_all(rf=rf)
        assert matrix.shape == (3003, 3003), "PWD matrix all failed!"


def test_similarity_matrix_all(ten_taxa_cMChain):
    for rf in [True, False]:
        matrix = ten_taxa_cMChain.similarity_matrix_all(rf=rf)
        assert matrix.shape == (3003, 3003), "Similarity matrix computation failed!"


def test_clustree_all(ten_taxa_cMChain):
    for rf in [True, False]:
        clustering = ten_taxa_cMChain.clustree_all(rf=rf)
        assert clustering.shape == (3003,), "Clustering all trees failed!"


def test_plot_clustree_all(ten_taxa_cMChain):
    for rf in [True, False]:
        ten_taxa_cMChain.plot_clustree_all(rf=rf)
        assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_2clustering_all{'_rf' if rf else ''}.png"),\
            "Plotting clustree all failed!"


def test_plot_clustree_all(ten_taxa_cMChain):
    for rf in [True, False]:
        ten_taxa_cMChain.plot_clustree_all(rf=rf, n_clus=1)
        assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_1clustering_all{'_rf' if rf else ''}.png"),\
            "Plotting clustree all failed!"


def test_tsne_all(ten_taxa_cMChain):
    for rf in [True, False]:
        ten_taxa_cMChain.tsne_all(rf=rf)
        assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/data/{ten_taxa_cMChain.name}_2D_tSNE_all{'_rf' if rf else ''}.npy"), \
            "2D tSNE all chains failed!"


def test_cen_for_each_chain(ten_taxa_cMChain):
    ten_taxa_cMChain.cen_for_each_chain()
    for chain in ten_taxa_cMChain:
        assert os.path.exists(f"{chain.working_dir}/{chain.name}_cen.tree"), "Failed to write centroid file!"


def test_compare_chain_summaries(ten_taxa_cMChain):
    try:
        os.remove(f"{ten_taxa_cMChain.working_dir}/data/{ten_taxa_cMChain.name}_cen_distances.log")
    except FileNotFoundError:
        pass
    ten_taxa_cMChain.compare_chain_summaries()

    # todo make list of TRUEs with os.path.exists instead of this one path
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/data/{ten_taxa_cMChain.name}_cen_distances.log"), "Compare chain summaries failed!"


def test_plot_chains_with_summaries(ten_taxa_cMChain):
    for rf in [True, False]:
        ten_taxa_cMChain.plot_chains_with_summaries(rf=rf)
        assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_cen_summaries{'_rf' if rf else ''}.png"), \
            "Plot chains with summaries failed!"


def test_compare_cutoff_treesets(ten_taxa_cMChain):
    for ess in [0, 200, 500]:
        ten_taxa_cMChain.compare_cutoff_treesets(i=0, j=1, ess=ess, beast_applauncher="/home/lars/.local/share/beast/bin/applauncher", _overwrite=True)
    assert True  # todo better test assertion for this one


def test_compare_cutoff_ess_choices(ten_taxa_cMChain):
    try:
        os.remove(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_0_1_[0]_cutoff_full_comparison.png")
    except FileNotFoundError:
        pass
    ten_taxa_cMChain.compare_cutoff_ess_choices(i=0, j=1)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_0_1_[0]_cutoff_full_comparison.png"), \
        "Compare Cutoff ESS choices plot failed!"


def test_compare_cutoff_ess_choices_all(ten_taxa_cMChain):
    try:
        os.remove(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_0_1_[0, 200, 500]_cutoff_full_comparison.png")
    except FileNotFoundError:
        pass
    ten_taxa_cMChain.compare_cutoff_ess_choices(i=0, j=1, ess_l=[0, 200, 500])
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_0_1_[0, 200, 500]_cutoff_full_comparison.png"), \
        "Compare Cutoff ESS all choices plot failed!"


def test_clade_set_comparison(ten_taxa_cMChain):
    burnin = 0.1
    try:
        os.remove(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_discrete_cc_0_1_burn-{burnin}.png")
    except FileNotFoundError:
        pass
    ten_taxa_cMChain.clade_set_comparison(i=0, j=1, burnin=burnin)
    assert os.path.exists(f"{ten_taxa_cMChain.working_dir}/plots/{ten_taxa_cMChain.name}_discrete_cc_0_1_burn-{burnin}.png"), "Discrete clade set comparator failed!"


def test_clade_set_comparison_noplot(ten_taxa_cMChain):
    v = ten_taxa_cMChain.clade_set_comparison(i=0, j=1, plot=False)
    assert v == 99.93756243756243, "Discrete clade set comparison without plot failed!"


#todo temporary
def test_cMChain_split_all_trees(ten_taxa_cMChain):
    ten_taxa_cMChain.split_all_trees(n_clus=3)
    # todo this test should first remove the files and then run the test
    # todo this is currently not implemented for rf option !!!
    assert True

def test_get_best_bic_cluster(ten_taxa_cMChain):
    best_bic = ten_taxa_cMChain[0].get_best_bic_cluster()
    assert best_bic == 5, "Get best bic Mchain function failed"


def test_get_best_silhouette_cluster(ten_taxa_cMChain):
    best_sil = ten_taxa_cMChain[0].get_best_silhouette_cluster()
    assert best_sil == 5, "Get best silhouette cluster evaluation failed"
