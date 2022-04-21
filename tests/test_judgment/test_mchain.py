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
    ),"Construction files failed"


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

