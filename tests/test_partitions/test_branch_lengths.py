from conpardist.branch_lengths import get_lengths


def test_get_conditional_branchlengths(thirty_taxa_tts):
    bl = get_lengths(thirty_taxa_tts[0])
    assert True

