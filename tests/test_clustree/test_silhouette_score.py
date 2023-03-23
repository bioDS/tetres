from tetres.clustree.silhouette_score import plot_silhouette, silhouette_score


def test_silhouette_score(ten_taxa_cMChain):
    ret = silhouette_score(matrix=ten_taxa_cMChain.pwd_matrix(0),
                           k=2, local_norm=False,
                           working_folder=ten_taxa_cMChain.working_dir,
                           random_shuffle=False
                           )
    assert ret == 0.17390331975243276, "Silhouette using sklearn failed!"


def test_silhouette_score_local(ten_taxa_cMChain):
    ret = silhouette_score(matrix=ten_taxa_cMChain.pwd_matrix(0),
                           k=2, local_norm=True,
                           working_folder=ten_taxa_cMChain.working_dir,
                           random_shuffle=False
                           )
    assert ret == 0.389222072562057, "Silhouette using local-norm failed!"


def test_plot_silhouette(ten_taxa_cMChain):
    plot_silhouette(matrix=ten_taxa_cMChain.pwd_matrix(0), local_norm=False,
                           working_folder=ten_taxa_cMChain.working_dir, _overwrite=True,
                           add_random=True, name=ten_taxa_cMChain[0].name)
    assert True


def test_plot_silhouette_local(ten_taxa_cMChain):
    plot_silhouette(matrix=ten_taxa_cMChain.pwd_matrix(0), local_norm=True,
                           working_folder=ten_taxa_cMChain.working_dir, _overwrite=True,
                           add_random=True, name=ten_taxa_cMChain[0].name)
    assert True
