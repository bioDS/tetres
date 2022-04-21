import pytest
import pandas as pd
from pathlib import Path
from treeoclock.judgment.mchain import _read_beast_logfile, MChain


def test_read_beast_logfile():
    log_file_path = f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log"
    assert type(_read_beast_logfile(log_file_path)) is pd.DataFrame, "Empty function failed!"


def test_MChain_construction_files():
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction files failed"


def test_MChain_construction_TimeTreeSets(thirty_taxa_tts, thirty_taxa_centroid_tts):
    assert MChain(
        trees=thirty_taxa_tts,
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=thirty_taxa_centroid_tts,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction TTS failed"


def test_MChain_construction_mixed1(thirty_taxa_tts):
    assert MChain(
        trees=thirty_taxa_tts,
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction TTS trees, centroid file failed"


def test_MChain_construciton_mixed2(thirty_taxa_centroid_tts):
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=thirty_taxa_centroid_tts,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction trees file, TTS centroid failed"


def test_MChain_construciton_summary_None():
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=None,
        working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
    ), "Construction with summary=None failed"


def test_MChain_construciton_workingdir_None():
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
        working_dir=None
    ), "Construction with working_dir=None failed"


def test_MChain_construciton_workingdir_summary_None():
    assert MChain(
        trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
        log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
        summary=None,
        working_dir=None
    ), "Construction with working_dir=None failed"


# testing construction exceptions
def test_MChain_wrong_trees():
    with pytest.raises(ValueError):
        MChain(
            trees=20,
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_wrong_logfile_path():
    with pytest.raises(FileNotFoundError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_wrong_logfile_type():
    with pytest.raises(ValueError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=20,
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_wrong_summary():
    with pytest.raises(ValueError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=20,
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_wrong_workingdir_path():
    with pytest.raises(NotADirectoryError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/work/"
        )


def test_MChain_wrong_workingdir_type():
    with pytest.raises(ValueError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_centroid.tree",
            working_dir=20
        )


def test_MChain_wrong_treemaps():
    with pytest.raises(ValueError):
        MChain(
            trees=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa.trees",
            log_file=f"{Path(__file__).parent.parent.absolute()}/data/30Taxa_beast2.log",
            summary=f"{Path(__file__).parent.parent.absolute()}/data/12Taxa.trees",
            working_dir=f"{Path(__file__).parent.parent.absolute()}/data/"
        )


def test_MChain_get_key_names(thirty_taxa_MChain):
    assert thirty_taxa_MChain.get_key_names() == ['posterior',
                                                  'likelihood',
                                                  'prior',
                                                  'treeLikelihood',
                                                  'TreeHeight',
                                                  'popSize',
                                                  'CoalescentConstant'], \
        "Get_key_name functions failed!"


def test_MChain_get_ess_tracerer(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key, "tracerer"))
    assert result == [1476.6549659664809,
                      1203.6766888338857,
                      1668.228049029279,
                      1203.6766888338857,
                      1758.881844170911,
                      1650.888969360209,
                      1677.680147437517], \
        "ESS tracerer value calculation with MChain failed!"


def test_MChain_get_ess_coda(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key, "coda"))
    assert result == [1447.80470846111,
                      1398.5896384868036,
                      1697.2455114765976,
                      1398.5896384868036,
                      1819.6601044141469,
                      1635.0719443144444,
                      1707.950568955924], \
        "ESS coda value calculation with MChain failed!"


def test_MChain_get_ess_arviz(thirty_taxa_MChain):
    result = []
    for ess_key in thirty_taxa_MChain.get_key_names():
        result.append(thirty_taxa_MChain.get_ess(ess_key, "arviz"))
    assert result == [1446.1461134383974,
                      1388.2577408057957,
                      1680.2897556479547,
                      1388.2577408057957,
                      1830.3012635184548,
                      1661.0004487633273,
                      1683.7924322905303], \
        "ESS arviz value calculation with MChain failed!"
