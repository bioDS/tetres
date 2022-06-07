from treeoclock.beast_helper.beast import run_beast
from pathlib import Path
import os


def test_run_beast():
    run_beast(xml_file=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/DS1.xml",
              beast="/home/lars/.local/share/beast/bin/beast")
    assert os.path.exists(f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/DS1.trees"),\
        "run_beast failed!"
