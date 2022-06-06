from treeoclock.beast_helper.beter import process_template, _init_R
from pathlib import Path
import pytest
import os


def test_process_template_wrongfiles_nowd():
    with pytest.raises(FileNotFoundError):
        process_template("DS1-template.xml", "DS1.xml", "config.toml")


def test_process_template_wrongfiles_wd():
    with pytest.raises(FileNotFoundError):
        process_template("DS1-template.xml", "DS1.xml", "config_1.toml",
                         working_dir=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper")


def test_process_template_wrongfiles_wd():
    with pytest.raises(FileNotFoundError):
        process_template("DS1-template_1.xml", "DS1.xml", "config.toml",
                         working_dir=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper")


def test_process_template_wrongwd():
    with pytest.raises(NotADirectoryError):
        process_template("DS1-template.xml", "DS1.xml", "config.toml",
                         working_dir=f"{Path(__file__).parent.parent.absolute()}/data/beast_help")


def test_process_template_wd():
    os.remove(f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/DS1.xml")
    process_template("DS1-template.xml", "DS1.xml", "config.toml",
                     working_dir=f"{Path(__file__).parent.parent.absolute()}/data/beast_helper")
    assert os.path.exists(f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/DS1.xml"), \
        "process_template failed with working_dir given!"


def test_process_template_no_wd():
    os.remove(f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/DS1.xml")
    working_dir = f"{Path(__file__).parent.parent.absolute()}/data/beast_helper"
    process_template(f"{working_dir}/DS1-template.xml", f"{working_dir}/DS1.xml",
                     f"{working_dir}/config.toml")
    assert os.path.exists(f"{Path(__file__).parent.parent.absolute()}/data/beast_helper/DS1.xml"), \
        "process_template failed without working_dir given!"
