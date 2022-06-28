from treeoclock.clustree.spectral_clustree import spectral_clustree_dm


def test_spectral_clustree(ten_taxa_cMChain):
    c = spectral_clustree_dm(ten_taxa_cMChain[0].pwd_matrix())
    print(c)
    assert True