from tetres.clustree.spectral_clustree import spectral_clustree_dm


def test_spectral_clustree(ten_taxa_multichain):
    c = spectral_clustree_dm(ten_taxa_multichain[0].pwd_matrix())
    print(c)
    assert True
