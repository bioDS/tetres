from rpy2.robjects.packages import importr, isinstalled
import os


def _init_R(func):
    if not isinstalled("beter"):
        raise ValueError("Package 'beter' is not installed! "
                         "Please do manually using 'devtools::install_github('biods/beter')'")
    return func


@_init_R
def process_template(xml_template, out_file, config, working_dir):
    beter = importr("beter")
    if working_dir is None:
        if not os.path.exists(xml_template):
            raise FileNotFoundError(f"XML template {xml_template} does not exist!")

        if not os.path.exists(config):
            raise FileNotFoundError(f"Config file {config} does not exist!")
        beter.process_template(f"{xml_template}", f"{out_file}", f"{config}")
    else:
        if not os.path.isdir(working_dir):
            raise NotADirectoryError(f"Working directory {working_dir} does not exist!")
        if not os.path.exists(f"{working_dir}/{xml_template}"):
            raise FileNotFoundError(f"XML template {xml_template} file does not exist!")
        if not os.path.exists(f"{working_dir}/{config}"):
            raise FileNotFoundError(f"Config file {config} does not exist!")

        beter.process_template(f"{working_dir}/{xml_template}",
                               f"{working_dir}/{out_file}",
                               f"{working_dir}/{config}")

