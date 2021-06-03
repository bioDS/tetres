from treeoclock.trees.time_trees import TimeTreeSet, TimeTree, findpath_distance
from treeoclock.summary.frechet_mean import frechet_mean


def test_frechet_mean(five_taxa_tts):
    fm = frechet_mean([five_taxa_tts[0], five_taxa_tts[2]])
    tfm = TimeTree('((4:2,(3:1,5:1):1):2,(1:3,2:3):1);')
    d = [findpath_distance(tfm, i) for i in [five_taxa_tts[0], five_taxa_tts[2], fm]]
    assert d == [1, 1, 0], 'Frechet Mean computation failed!'
