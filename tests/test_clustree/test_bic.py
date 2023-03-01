from treeoclock.clustree.bic import BIC, plot_bic
import os


def test_BIC_1_false(ten_taxa_cMChain):

    bic, mnd_cluster = BIC(treeset=ten_taxa_cMChain[0].trees,
                           matrix=ten_taxa_cMChain.pwd_matrix(0),
                           k=1, local=False,
                           working_folder=ten_taxa_cMChain.working_dir,
                           add_random=False,
                           matrix_is_similarity=False
                           )
    assert int(bic) == 395, "BIC for 1 cluster failed"


def test_BIC_1_true(ten_taxa_cMChain):
    bic, mnd_cluster = BIC(treeset=ten_taxa_cMChain[0].trees,
                           matrix=ten_taxa_cMChain.pwd_matrix(0),
                           k=2, local_norm=False,
                           working_folder=ten_taxa_cMChain.working_dir,
                           add_random=False, _overwrite=True)
    assert int(bic) == 398, "BIC for 2 clusters failed"


def test_plot_bic(ten_taxa_cMChain):
    plot_bic(treeset=ten_taxa_cMChain[0].trees,
             matrix=ten_taxa_cMChain.pwd_matrix(0),
             working_folder=ten_taxa_cMChain.working_dir,
             name=ten_taxa_cMChain[0].name,
             _overwrite=True
             )

def test_plot_bic_reading(ten_taxa_cMChain):
    try:
        os.remove(f"{ten_taxa_cMChain.working_dir}/plots/BIC_10TaxaFix_0.png")
    except FileNotFoundError:
        pass
    plot_bic(treeset=ten_taxa_cMChain[0].trees,
             matrix=ten_taxa_cMChain.pwd_matrix(0),
             working_folder=ten_taxa_cMChain.working_dir,
             name=ten_taxa_cMChain[0].name,
             _overwrite=False
             )


def test_plot_bic_local_scale(ten_taxa_cMChain):
    try:
        os.remove(f"{ten_taxa_cMChain.working_dir}/plots/BIC_local_10TaxaFix_0.png")
    except FileNotFoundError:
        pass
    plot_bic(treeset=ten_taxa_cMChain[0].trees,
             matrix=ten_taxa_cMChain.pwd_matrix(0),
             working_folder=ten_taxa_cMChain.working_dir,
             name=ten_taxa_cMChain[0].name,
             _overwrite=False,
             local_norm=True
             )


def test_plot_bic_add_random(ten_taxa_cMChain):
    try:
        os.remove(f"{ten_taxa_cMChain.working_dir}/plots/BIC_random_10TaxaFix_0.png")
    except FileNotFoundError:
        pass
    plot_bic(treeset=ten_taxa_cMChain[0].trees,
             matrix=ten_taxa_cMChain.pwd_matrix(0),
             working_folder=ten_taxa_cMChain.working_dir,
             name=ten_taxa_cMChain[0].name,
             _overwrite=False,
             add_random=True
             )
