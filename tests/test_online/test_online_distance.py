
from treeoclock.online.online_distance import online_distance_neighbourhood


def test_online_distance_neighbourhood(twelve_taxa_tts):
    assert online_distance_neighbourhood(twelve_taxa_tts[0], twelve_taxa_tts[1]) == 0,\
        "Online distance neighbourhood failed!"
