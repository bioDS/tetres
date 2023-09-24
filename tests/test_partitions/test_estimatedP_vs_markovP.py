from tetres.partitions.estimatedP_vs_markovP import *

def test_ccd_cpd_vs_markov(twenty_taxa_tts):
    ud, rd = ccd_cpd_vs_markov(twenty_taxa_tts)
    assert (len(ud), len(rd)) == (1219, 2499), "CCD, CPD data went wrong."
