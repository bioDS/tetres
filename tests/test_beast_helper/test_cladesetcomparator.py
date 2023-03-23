from tetres.beast_helper.cladesetcomparator import run_cladesetcomparator
from pathlib import Path
import os

def test_run_beast():
    run_cladesetcomparator(tree_file1=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/run1/DS1.trees",
                            tree_file2=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/run2/DS1.trees",
                            out_file_png=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/r1_r2_cc.png",
                            beast_applauncher="/home/lars/.local/share/beast/bin/applauncher")

    assert os.path.exists(f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/r1_r2_cc.png"),\
        "run_cladesetcomparator failed!"
