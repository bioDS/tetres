import pytest

@pytest.fixture
def five_taxa_newick_list():
    # List of two 5 taxa trees with a distance of 1
    return ["((1:3,5:3):1,(4:2,(3:1,2:1):1):2);", "((1:2,5:2):3,(4:3,(3:1,2:1):2):2);"]
