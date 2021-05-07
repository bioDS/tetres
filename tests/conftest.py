import pytest

@pytest.fixture
def five_taxa_newick_list():
    # list of pairwise distances is [0, 1, 2, 1, 0, 3, 2, 3, 0]
    return ["((1:3,5:3):1,(4:2,(3:1,2:1):1):2);", "((1:2,5:2):3,(4:3,(3:1,2:1):2):2);",
            "((1:2,5:2):3,(3:3,(4:1,2:1):2):2);"]
