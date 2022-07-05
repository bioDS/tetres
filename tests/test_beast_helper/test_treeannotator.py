from treeoclock.beast_helper.treeannotator import run_treeannotator
from pathlib import Path
import os


def test_treeannotator():
    run_treeannotator(tree_file=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/run1/DS1.trees",
                      outfile=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/run1/mcc.tree",
                      burnin=5,
                      beast_treeannotator="/home/lars/.local/share/beast/bin/treeannotator")
    assert os.path.exists(f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/run1/mcc.tree"),\
        "Treeannotator helper failed!"
