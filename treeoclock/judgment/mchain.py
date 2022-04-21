import pandas as pd
import os

from treeoclock.judgment import ess
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
            self.chain_length = int(list(self.log_data["Sample"])[-1])
            self.sampling_interval = int(list(self.log_data["Sample"])[1])
        else:
            if type(log_file) is str:
                raise FileNotFoundError(log_file)
            else:
                raise ValueError(log_file)

        if summary is not None:
            # Setting the summary tree
            if type(summary) is TimeTreeSet:
                self.summary = summary
            elif type(summary) is str:
                self.summary = TimeTreeSet(summary)
            else:
                raise ValueError(summary)

        if working_dir is not None:
            # Setting the Working directory
            if os.path.isdir(working_dir):
                self.working_dir = working_dir
            else:
                if type(working_dir) is str:
                    raise NotADirectoryError(working_dir)
                else:
                    raise ValueError(working_dir)

        if summary is not None:
            # Check mapping of summary tree and tree set
            if not self.summary.map == self.trees.map:
                try:
                    self.summary.change_mapping(self.trees.map)
                except ValueError as error:
                    raise ValueError(f"{error}\n"
                                     f"The given summary tree and tree set do not fit! "
                                     f"\n(Construction of class MChain failed!)")

    def get_key_names(self):
        return list(self.log_data.columns)[1:]

    def get_ess(self, ess_key, ess_method):
        if type(ess_method) is str:
            if not hasattr(ess, f"{ess_method}_ess"):
                raise ValueError(f"The given ESS method {ess_method} does not exist!")
        else:
            raise ValueError(ess_method)

        if ess_key in list(self.log_data.columns):
            return getattr(ess, f"{ess_method}_ess")(data_list=self.log_data[ess_key], chain_length=self.chain_length, sampling_interval=self.sampling_interval)
        else:
            raise ValueError("Not (yet) implemented!")

    def get_ess_partial(self, ess_key, ess_method, lower_i=0, upper_i=1):
        # todo returns the ess value of a part of the MChain object, by default it returns the same es get_ess
        return 0

    def get_ess_trace_plot(self):
        # todo needs a list of ess_keys and ess_method
        #  interval size defaulting to 1 so all samples
        return 0

    # todo change the name of this plot!!!!
    def get_burnin_comparison(self):
        # todo gets a focal tree defaulting to summary tree
        #  if that is None, either pick one of the trees or compute something with options here
        return 0


def _read_beast_logfile(logfile_path):
    data = pd.read_csv(logfile_path, header=0, sep=r"\s+", comment="#", error_bad_lines=True)
    return data
