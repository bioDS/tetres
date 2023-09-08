from tetres.judgement.burnin_detection import burn_detector
import pytest


def test_burn_detector(ten_taxa_multichain):
    burn_detector(ten_taxa_multichain, 0, 1)
    assert True


def test_burn_detector_valueerror(ten_taxa_multichain):
    with pytest.raises(ValueError):
        burn_detector(ten_taxa_multichain, 0, 1, traces = ["posterior","blah"])
