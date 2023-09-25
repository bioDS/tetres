from tetres.partitions.estimatedP_vs_markovP import *

def test_ccd_cpd_vs_markov(twenty_taxa_tts):
    ud, rd = ccd_cpd_vs_markov(twenty_taxa_tts)
    assert (len(ud), len(rd)) == (1219, 2499), "CCD, CPD data went wrong."


def test_ccd_cpd_vs_markov_norank(twenty_taxa_tts):
    ud = ccd_cpd_vs_markov(twenty_taxa_tts, include_rank=False)
    assert len(ud) == 1219, "CCD, CPD data went wrong."
