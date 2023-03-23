from tetres.visualize.tsne import tsne_coords_from_mchain


def test_tsne_from_mchain_2d(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[0], dim=2)
    assert len(coords) == len(ten_taxa_cMChain[0].trees), "2D tsne failed!"


def test_tsne_from_mchain_3d(ten_taxa_cMChain):
    coords = tsne_coords_from_mchain(ten_taxa_cMChain[0], dim=3)
    assert len(coords) == len(ten_taxa_cMChain[0].trees), "3D tsne failed!"
