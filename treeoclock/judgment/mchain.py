import pandas as pd
import os

from treeoclock.trees.time_trees import TimeTreeSet


class MChain:
    def __init__(self, trees, log_file, summary, working_dir):

        # Setting trees
        if type(trees) is TimeTreeSet:
            self.trees = trees
        elif type(trees) is str:
            self.trees = TimeTreeSet(trees)
        else:
            raise ValueError(trees)

        # Reading the log_file
        if os.path.exists(log_file):
            self.log_data = _read_beast_logfile(log_file)
        else:
            if type(log_file) is str:
                raise FileNotFoundError(log_file)
            else:
                raise ValueError(log_file)

        # Setting the summary tree
        if type(summary) is TimeTreeSet:
            self.summary = summary
        elif type(summary) is str:
            self.summary = TimeTreeSet(summary)
        else:
            raise ValueError(summary)

        # Setting the Working directory
        if os.path.isdir(working_dir):
            self.working_dir = working_dir
        else:
            if type(working_dir) is str:
                raise NotADirectoryError(working_dir)
            else:
                raise ValueError(working_dir)
        if not self.summary.map == self.trees.map:
            try:
                self.summary.change_mapping(self.trees.map)
            except ValueError as error:
                raise ValueError(f"{error}\n"
                                 f"The given summary tree and tree set do not fit! "
                                 f"\n(Construction of class MChain failed!)")


def _read_beast_logfile(logfile_path):
    data = pd.read_csv(logfile_path, header=0, sep=r"\s+", comment="#", error_bad_lines=True)
    return data
