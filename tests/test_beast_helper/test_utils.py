import os

from tetres.beast_helper.utils import change_taxon_map


def test_change_taxon_map(ten_taxa_multichain):
    in_nexus = os.path.join(ten_taxa_multichain.working_dir, ten_taxa_multichain.tree_files[0])
    out_nexus = os.path.join(ten_taxa_multichain.working_dir, f"changed_mapping_{ten_taxa_multichain.tree_files[0]}")
    new_map = { 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10}
    change_taxon_map(input_nexus=in_nexus, output_nexus=out_nexus, new_map=new_map)
    assert True


def test_change_taxon_map_mcc(ten_taxa_multichain):
    in_nexus = os.path.join(ten_taxa_multichain.working_dir, f"chain0_1-mcc.tree")
    out_nexus = os.path.join(ten_taxa_multichain.working_dir, f"changed_mapping_mcc.tree")
    new_map = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10}
    change_taxon_map(input_nexus=in_nexus, output_nexus=out_nexus, new_map=new_map)
    assert True
