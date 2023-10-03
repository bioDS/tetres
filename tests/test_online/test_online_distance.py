
from tetres.online.online_distance import online_distance_neighbourhood, new_function


def test_online_distance_neighbourhood(twelve_taxa_tts):
    d = online_distance_neighbourhood(twelve_taxa_tts[0], twelve_taxa_tts[1])
    assert True


def test_new_c_stuff(twelve_taxa_tts):
    res = new_function(twelve_taxa_tts)
    assert res == 1, "New Stuff failed"
