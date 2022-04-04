from treeoclock.trees.time_trees import TimeTree, findpath_distance
from treeoclock.summary.frechet_mean import frechet_mean, frechet_mean_sort
from treeoclock.summary.compute_sos import compute_sos

import random

def test_frechet_mean(five_taxa_tts):
    fm = frechet_mean([five_taxa_tts[0], five_taxa_tts[2]])
    tfm = TimeTree('((4:2,(3:1,5:1):1):2,(1:3,2:3):1);')
    d = [findpath_distance(tfm, i) for i in [five_taxa_tts[0], five_taxa_tts[2], fm]]
    assert d == [1, 1, 0], 'Frechet Mean computation failed!'


def test_frechet_mean_data(twelve_taxa_tts):
    state = random.getstate()  # get the random seed state
    random.seed(10)  # Fixing the seed to get the same result
    fm = frechet_mean(twelve_taxa_tts)
    sos = compute_sos(fm, twelve_taxa_tts)
    random.setstate(state)  # reset the random seed to previous state
    assert sos == 351747, "Frechet Mean failed on data set!"


def test_frechet_mean_sort(twelve_taxa_tts):
    fm = frechet_mean_sort(twelve_taxa_tts)
    sos = compute_sos(fm, twelve_taxa_tts)
    assert sos == 289661, "Frechet mean sort failed!"