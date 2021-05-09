import pytest

from testfixtures import TempDirectory


@pytest.fixture()
def dir():
    with TempDirectory() as dir:
        yield dir


@pytest.fixture
def five_taxa_newick_list():
    # list of pairwise distances is [0, 1, 2, 1, 0, 3, 2, 3, 0]
    return ["((1:3,5:3):1,(4:2,(3:1,2:1):1):2);", "((1:2,5:2):3,(4:3,(3:1,2:1):2):2);",
            "((1:2,5:2):3,(3:3,(4:1,2:1):2):2);"]


@pytest.fixture
def five_taxa_newick_list_distances():
    return [0, 1, 2, 1, 0, 3, 2, 3, 0]


@pytest.fixture
def five_taxa_nexus_string():
    return ("#NEXUS\n"
            "BEGIN TAXA;\n"
            "\tDIMENSIONS NTAX = 5;\n"
            "\tTAXLABELS\n"
            "\t\tt1\n"
            "\t\tt5\n"
            "\t\tt4\n"
            "\t\tt3\n"
            "\t\tt2\n"
            "\t;\n"
            "END;\n"
            "BEGIN TREES;\n"
            "\tTRANSLATE\n"
            "\t1	t1,\n"
            "\t2	t5,\n"
            "\t3	t4,\n"
            "\t4	t3,\n"
            "\t5	t2\n"
            "\t;\n"
            "\tTREE * t1 = [&R] ((1:3,2:3):1,(3:2,(4:1,5:1):1):2);\n"
            "\tTREE * t2 = [&R] ((1:2,2:2):3,(3:3,(4:1,5:1):2):2);\n"
            "\tTREE * t3 = [&R] ((1:2,2:2):3,(4:3,(3:1,5:1):2):2);\n"
            "END;\n").encode()
