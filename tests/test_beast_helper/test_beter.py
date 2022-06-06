from treeoclock.beast_helper.beter import process_template, _init_R
from pathlib import Path


def test_process_template():
    process_template("DS1-template.xml", "DS1.xml", "config.toml",
                     working_dir=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper")
    assert True

# todo tests with wrong files given to trigger exceptions