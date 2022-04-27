import sys
import os
import pandas as pd
import numpy as np


from treeoclock.judgment import ess, _ess_plots
from treeoclock.trees.time_trees import TimeTreeSet
from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary.frechet_mean import frechet_mean


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
            self.log_sampling_interval = int(list(self.log_data["Sample"])[1])
            # assuming that the first iteration 0 tree and the last have been logged in the logfile
            self.tree_sampling_interval = self.chain_length / (len(self.trees) - 1)
        else:
            if type(log_file) is str:
                raise FileNotFoundError(log_file)
            else:
                raise ValueError(log_file)

        # todo future work for this part!
        # if self.log_data.shape[0] != len(self.trees):
        #     raise ValueError("Different Size of Tree and LogFile!")

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

    def get_ess(self, ess_key="posterior", ess_method="arviz", **kwargs):
        if type(ess_method) is str:
            if not hasattr(ess, f"{ess_method}_ess"):
                raise ValueError(f"The given ESS method {ess_method} does not exist!")
        else:
            raise ValueError(ess_method)
        lower_i = 0
        if "lower_i" in kwargs:
            lower_i = kwargs["lower_i"]
        upper_i = self.log_data.shape[0]
        if "upper_i" in kwargs:
            upper_i = kwargs["upper_i"]
        if "lower_i" in kwargs or "upper_i" in kwargs:
            if type(lower_i) is not int or type(upper_i) is not int:
                raise ValueError("Wrong type for upper or lower index!")
            if upper_i > self.log_data.shape[0]:
                raise IndexError(f"{upper_i} out of range!")
            if lower_i > upper_i or lower_i < 0 or upper_i < 0:
                raise ValueError("Something went wrong with the given upper and lower index!")

        if ess_key in list(self.log_data.columns):
            chain_length = 1
            if upper_i != lower_i:
                chain_length = int(self.log_data["Sample"][upper_i])-int(self.log_data["Sample"][lower_i])
            cur_sampling_interval = int(chain_length / (self.log_data[ess_key][lower_i:(upper_i+1)].shape[0] - 1 - self.log_data[ess_key][lower_i:(upper_i+1)].isna().sum()))
            return getattr(ess, f"{ess_method}_ess")(data_list=self.log_data[ess_key][lower_i:(upper_i+1)].dropna(),
                                                     chain_length=chain_length,
                                                     sampling_interval=cur_sampling_interval)
        else:
            raise ValueError("Not (yet) implemented!")

    # todo ideally this should be accessible with the get_ess() funciton and ess_key="pseudo"
    def get_pseudo_ess(self, **kwargs):
        # todo all the testing for upper_i and lower_i but ideally this will be part of the other funciton so this is temporary
        # todo should return the same as RWTY package but also possibly computed with the RNNI distance measure
        upper_i = len(self.trees)
        if "upper_i" in kwargs:
            upper_i = kwargs["upper_i"]
        chain_length = 1
        if upper_i > 0:
            chain_length = (self.chain_length / (len(self.trees) - 1)) * (upper_i - 1)
        return ess.pseudo_ess(tree_set=self.trees[0:upper_i], chain_length=chain_length,
                              sampling_interval=self.chain_length / (len(self.trees) - 1))

    def compute_rnni_variance_log(self, focal_tree_type="tree", add=True):
        new_log_list = []
        for i in range(0, len(self.trees)):
            # todo this type checking is not efficient, should be outside of the loop!
            if focal_tree_type == "FM":
                focal_tree = frechet_mean(self.trees[0:i+1])
            elif focal_tree_type == "tree":
                focal_tree = self.trees[i]
            elif focal_tree_type == "centroid":
                sys.exit("Not yet implemented")
                focal_tree = frechet_mean(self.trees[0:i + 1])
            else:
                raise ValueError(f"Unknown type {focal_tree_type} given!")
            new_log_list.append(compute_sos_mt(focal_tree, self.trees[0:i+1], n_cores=None)/(i+1))

        # adding the computed log value to the log dataframe

        if add:
            if len(new_log_list) == self.log_data.shape[0]:
                self.log_data[f"Var_{focal_tree_type}"] = new_log_list
            elif len(new_log_list) < self.log_data.shape[0]:
                if self.tree_sampling_interval % self.log_sampling_interval == 0:
                    sampling_diff = int(self.tree_sampling_interval / self.log_sampling_interval)
                    self.log_data[f"Var_{focal_tree_type}"] = np.nan
                    self.log_data[f"Var_{focal_tree_type}"][::sampling_diff] = new_log_list
                else:
                    sys.exit("Feature not implemented! Trees logged is not a multiple of the log sampling interval!")
            elif len(new_log_list) > self.log_data.shape[0]:
                sys.exit("Feature not yet implemented! "
                         "The logfile contains less data than the tree file!")
        return new_log_list

    def get_ess_trace_plot(self, ess_key="all", ess_method="arviz", kind="cummulative"):
        # todo add kind window ?
        #  add a interval size for the plot
        # todo all the checks for the variables given
        # todo ess method can be list of all methods and then it will create multiplot thing with distribution on the diagonal

        if ess_key == "all":
            ess_key = self.get_key_names()
        elif type(ess_key) is str and ess_key in self.log_data:
            ess_key = [ess_key]
        elif not (type(ess_key) is list and all(item in self.get_key_names() for item in ess_key)):
            raise ValueError("ess_key wrong type")

        data = []
        for i in range(5, self.log_data.shape[0]):  # Starting at sample 5 as it is not useful to look at less samples
            data.extend([[key, self.get_ess(ess_key=key, ess_method=ess_method, upper_i=i), i] for key in ess_key])
        data = pd.DataFrame(data, columns=["Ess_key", "Ess_value", "Upper_i"])
        _ess_plots._ess_trace_plot(data)
        return 0

    # todo change the name of this plot!!!!
    def get_burnin_comparison(self):
        # todo gets a focal tree defaulting to summary tree
        #  if that is None, either pick one of the trees or compute something with options here
        return 0


def _read_beast_logfile(logfile_path):
    data = pd.read_csv(logfile_path, header=0, sep=r"\s+", comment="#", error_bad_lines=True)
    return data
