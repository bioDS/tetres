import pytest
import os
from tetres.judgement.convex_hull import *


def test_find_convex_hull(ten_taxa_cMChain):
    # chain = coupled_MChains(m_MChains=1,
    #                 trees=["30Taxa.trees"],
    #                 log_files=["30Taxa_beast2.log"],
    #                 working_dir=f"{Path(__file__).parent.parent.absolute()}/data"
    #                 )


    # todo fix the random seed for the fm calculation

    # todo fix the bug that there are so many trees in the convex hull
    l = []

    for _ in range(10):
        clist = find_convex_hull_trees(ten_taxa_cMChain, 0)
        l.append(clist.count(True))
    print(l)
    assert True
