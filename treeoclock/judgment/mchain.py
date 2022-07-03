import sys
import os
import pandas as pd
import numpy as np
import itertools
import linecache
import re

from treeoclock.judgment import ess, _plots
from treeoclock.judgment import _geweke_diag as gwd
from treeoclock.trees.time_trees import TimeTreeSet, TimeTree
from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary.frechet_mean import frechet_mean
from treeoclock import enums

from treeoclock.judgment._pairwise_distance_matrix import calc_pw_distances, calc_pw_distances_two_sets
from treeoclock.judgment import _gelman_rubin_diag as grd
from treeoclock.judgment import _cladesetcomparator as csc
from treeoclock.clustree.spectral_clustree import _spectral_clustree
from treeoclock.judgment._plotting import all_chains_spectral_clustree


# todo testing required!
class coupled_MChains():
    def __init__(self, m_MChains, trees, log_files, working_dir, name="cMC"):
        # todo currently still missing the summary parameter of MChain, but its not really used at this point
        self.m_MChains = m_MChains
        self.name = name

        if type(trees) != type(log_files):
            raise ValueError("trees and logfiles should be given with the same data type!")

        if not os.path.isdir(working_dir):
            raise FileNotFoundError(f"Given Working directory does not exist! |{working_dir}|")
        self.working_dir = working_dir
        for d in ["data", "plots"]:
            if not os.path.isdir(f"{self.working_dir}/{d}"):
                os.mkdir(f"{self.working_dir}/{d}")

        self.MChain_list = []
        if type(trees) is list:
            if len(trees) != self.m_MChains or len(log_files) != self.m_MChains:
                raise ValueError("m_chains and number of trees do not match!")
            for i in range(self.m_MChains):
                self.tree_files = trees.copy()
                self.MChain_list.append(MChain(trees=trees[i], log_file=log_files[i],
                                               working_dir=working_dir, name=f"{self.name}_{i}"))
        elif type(trees) is str:
            if not os.path.exists(f"{self.working_dir}/{trees}"):
                raise FileNotFoundError(f"Given trees file {self.working_dir}/{trees} does not exist!")
            if not os.path.exists(f"{self.working_dir}/{log_files}"):
                raise FileNotFoundError(f"Given trees file {self.working_dir}/{log_files} does not exist!")
            self.MChain_list.append(MChain(trees=trees, log_file=log_files, working_dir=working_dir, name=f"{self.name}_{0}"))
            for i in range(1, self.m_MChains):
                self.MChain_list.append(
                    MChain(trees=f"chain{i}{trees}", log_file=f"chain{i}{log_files}",
                           working_dir=working_dir, name=f"{self.name}_{i}"))
        else:
            raise ValueError("Unrecognized argument types of trees and log_files!")

    def pwd_matrix(self, index1, index2=None, csv: bool = False):
        if type(index1) is not int:
            raise ValueError("Unrecognized index type!")
        if index1 >= self.m_MChains:
            raise IndexError("Given Index out of range!")

        if index2 is None:
            return self.MChain_list[index1].pwd_matrix(index=index1, name=self.name, csv=csv)
        else:
            if type(index2) is not int:
                raise ValueError("Unrecognized index type!")
            if index1 == index2:
                raise ValueError("Indeces are equal, only call with one index for pairwise distances of one set!")
            if index2 >= self.m_MChains:
                raise IndexError("Given Index out of range!")
            if index2 < index1:
                index1, index2 = index2, index1
            if not os.path.exists(f"{self.working_dir}/data/{self.name}_{index1}_{index2}.{'csv.gz' if csv else 'npy'}"):
                dm = calc_pw_distances_two_sets(self.MChain_list[index1].trees, self.MChain_list[index2].trees)
                if csv:
                    np.savetxt(fname=f"{self.working_dir}/data/{self.name}_{index1}_{index2}.csv.gz", X=dm, delimiter=',',
                               fmt='%i')
                else:
                    np.save(file=f"{self.working_dir}/data/{self.name}_{index1}_{index2}.npy", arr=dm)
                return dm
            else:
                if csv:
                    return np.genfromtxt(fname=f"{self.working_dir}/data/{self.name}_{index1}_{index2}.csv.gz",
                                         delimiter=',', dtype=int)
                else:
                    return np.load(f"{self.working_dir}/data/{self.name}_{index1}_{index2}.npy")

    # todo testing required
    def pwd_matrix_all(self):
        if not os.path.exists(f"{self.working_dir}/data/{self.name}_all.npy"):
            combined_matrix = np.array([])
            for i in range(self.m_MChains):
                cur_row = np.array([])
                for j in range(self.m_MChains):
                    # print(i, j)
                    if i < j:
                        cur_row = np.concatenate((cur_row, self.pwd_matrix(i, j)),
                                                 axis=1) if cur_row.size else self.pwd_matrix(i, j)
                    elif i > j:
                        cur_row = np.concatenate((cur_row, np.zeros(self.pwd_matrix(j, i).shape)),
                                                 axis=1) if cur_row.size else np.zeros(self.pwd_matrix(j, i).shape)
                    elif i == j:
                        cur_row = np.concatenate((cur_row, self.pwd_matrix(i)),
                                                 axis=1) if cur_row.size else self.pwd_matrix(i)
                    # print(cur_row.shape)
                combined_matrix = np.concatenate((combined_matrix, cur_row)) if combined_matrix.size else cur_row
            np.save(file=f"{self.working_dir}/data/{self.name}_all.npy", arr=combined_matrix)
            return combined_matrix
        else:
            return np.load(f"{self.working_dir}/data/{self.name}_all.npy")

    def similarity_matrix_all(self, beta=1):
        if not os.path.exists(f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity_all.npy"):
            matrix = self.pwd_matrix_all()
            similarity = np.exp(-beta * matrix / matrix.std())
            np.save(file=f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity_all.npy",
                    arr=similarity)
            return similarity
        else:
            return np.load(f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity_all.npy")

    def clustree_all(self, n_clus=2, beta=1):
        if not os.path.exists(
                f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all.npy"):
            similarity = self.similarity_matrix_all(beta=beta)
            similarity = similarity + similarity.transpose()
            clustering = _spectral_clustree(similarity, n_clus=n_clus)
            np.save(
                file=f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all.npy",
                arr=clustering)
            return clustering
        else:
            return np.load(
                f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all.npy")

    def spectral_cluster_all(self, n_clus=2, beta=1):
        all_chains_spectral_clustree(self, beta=beta, n_clus=n_clus)

    def split_all_trees(self, n_clus, beta=1, burn_in=5):
        clustering = self.clustree_all(n_clus=n_clus, beta=beta)
        # todo split log and tree files into new files to assess them as individual analyses

        try:
            os.mkdir(f"{self.working_dir}/{n_clus}_cluster")
        except FileExistsError:
            # raise FileExistsError("Cluster directory already exists! No overwriting implemented!")
            pass

        # todo this split needs to be based on the runs as well, i.e. not global but for each chain the split needs to happen!

        sample_dict = {}
        for c in range(n_clus):
            try:
                f = open(f"{self.working_dir}/{n_clus}_cluster/clus_{c}.log", "x")
                f.write("\t".join(v for v in list(self[0].log_data.keys())))
                f.write("\n")
                f.close()
                sample_dict[c] = 0
            except FileExistsError:
                raise FileExistsError("Log files from clustering already exist, manual deletion for recomputation "
                                      "required!")
            # tree files are initialized with the map and properties of chain at postion 0
            try:
                f = open(f"{self.working_dir}/{n_clus}_cluster/clus_{c}.trees", "x")
                f.write(f"#NEXUS\n\nBegin taxa;\n\tDimensions ntax={self[0].trees[0].ctree.num_leaves};\n\t\tTaxlabels\n")
                for taxa in range(1, self[0].trees[0].ctree.num_leaves+1):
                    f.write(f"\t\t\t{taxa}\n")
                f.write("\t\t\t;\nEnd;\nBegin trees;\n\tTranslate\n")
                for taxa in range(1, self[0].trees[0].ctree.num_leaves):
                    f.write(f"\t\t\t{taxa} {self[0].trees.map[taxa]},\n")
                f.write(f"\t\t\t{self[0].trees[0].ctree.num_leaves} {self[0].trees.map[self[0].trees[0].ctree.num_leaves]}\n")
                f.write(";\n")
                f.close()

            except FileExistsError:
                raise FileExistsError("Tree files from clustering already exist, manual deletion for recomputation "
                                      "required!")

        for chain in range(self.m_MChains):
            if chain != 0:
                if self[chain].trees.map != self[0].trees.map:
                    raise ValueError("Problem of unequal maps between the tree sets!")
            for index, row in self[chain].log_data.iterrows():
                if index > int((self[chain].log_data.shape[0] / 100) * burn_in):
                    # writing the log file
                    cur_file = open(f"{self.working_dir}/{n_clus}_cluster/clus_{clustering[index]}.log", "a")
                    cur_value_list = row.values
                    cur_value_list[0] = sample_dict[clustering[index]]
                    cur_file.write("\t".join([str(v) for v in cur_value_list]))
                    cur_file.write("\n")
                    cur_file.close()

                    # writing the tree file
                    cur_file = open(f"{self.working_dir}/{n_clus}_cluster/clus_{clustering[index]}.trees", "a")

                    re_tree = re.compile("\t?tree .*=? (.*$)", flags=re.I | re.MULTILINE)

                    offset = 10 + (2 * self[0].trees[0].ctree.num_leaves)
                    file = f"{self.working_dir}/{self.tree_files[chain]}"
                    line = index + 1 + offset
                    cur_tree = linecache.getline(file, line)
                    cur_tree = f'{re.split(re_tree, cur_tree)[1][:re.split(re_tree, cur_tree)[1].rfind(")") + 1]};'

                    if not TimeTree(cur_tree).fp_distance(self[chain].trees[index]) == 0:
                        raise ValueError("Wrong offset!")

                    cur_file.write(f"tree STATE_{sample_dict[clustering[index]]} = {cur_tree}\n")

                    cur_file.close()
                    # increasing the sample for the respective clustering
                    sample_dict[clustering[index]] += 1
        for c in range(n_clus):
            cur_file = open(f"{self.working_dir}/{n_clus}_cluster/clus_{c}.trees", "a")
            cur_file.write("End;")
            cur_file.close()
        return 0

    def gelman_rubin_like_diagnostic_plot(self, samples: int = 100):
        if len(self) < 2:
            raise ValueError("Gelman Rubin only possible with multiple Chains!")
        return grd.gelman_rubin_distance_diagnostic_plot(self, samples=samples)

    def gelman_rubin_trace_plot(self, i, j):
        if len(self) < 2:
            raise TypeError("Gelman Rubin only possible with multiple Chains!")
        return grd.gelman_rubin_trace_plot(self, i, j)

    def gelman_rubin_trace_ess_plot(self, i, j):
        if len(self) < 2:
            raise TypeError("Gelman Rubin only possible with multiple Chains!")
        return grd.gr_trace_ess(self, i, j)

    def __getitem__(self, index):
        if isinstance(index, int):
            if index >= len(self.MChain_list):
                raise IndexError("Index out of range!")
            return self.MChain_list[index]
        else:
            raise TypeError("Given Type for getitem not implemented!")

    def __len__(self):
        return len(self.MChain_list)

    def cladesetcomparator(self, beast_applauncher, burnin=10):
        if not self.tree_files:
            raise ValueError("Missing tree_files list!")
        csc._cladesetcomp(self, beast_applauncher, burnin=burnin)

    def ess_stripplot(self, ess_method="tracerer"):
        ess.ess_stripplot(self, ess_method)


class MChain:
    def __init__(self, working_dir, trees, log_file, summary=None, name: str = "MC"):
        if type(name) is str:
            self.name = name
        else:
            raise ValueError("Unrecognized type for name!")
        # Setting the Working directory
        if type(working_dir) is str:
            if os.path.isdir(working_dir):
                self.working_dir = working_dir
            else:
                if type(working_dir) is str:
                    raise NotADirectoryError(working_dir)
                else:
                    raise ValueError(working_dir)
        else:
            raise ValueError("Working directory type not recognized!")

        for d in ["data", "plots"]:
            if not os.path.isdir(f"{self.working_dir}/{d}"):
                os.mkdir(f"{self.working_dir}/{d}")

        # Setting trees
        if type(trees) is TimeTreeSet:
            self.trees = trees
        elif type(trees) is str:
            self.trees = TimeTreeSet(f"{self.working_dir}/{trees}")
        else:
            raise ValueError(trees)

        # Reading the log_file
        if os.path.exists(f"{self.working_dir}/{log_file}"):
            self.log_data = pd.read_csv(f"{self.working_dir}/{log_file}", header=0, sep=r"\s+", comment="#")
            self.chain_length = int(list(self.log_data["Sample"])[-1])
            self.log_sampling_interval = int(list(self.log_data["Sample"])[1])
            # assuming that the first iteration 0 tree and the last have been logged in the logfile
            self.tree_sampling_interval = int(self.chain_length / (len(self.trees) - 1))
        # todo the logfile should be optional as it is possible to compute all new measures without the log file
        #  it is only relevant if we want to compare them or do something with those values
        #  maybe add the option logfile None which also throws errors when functions are called that might use the logfile
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
                if os.path.exists(f"{self.working_dir}/{summary}"):
                    self.summary = TimeTreeSet(f"{self.working_dir}/{summary}")
                else:
                    raise FileNotFoundError("Given summary tree not found!")
            else:
                raise ValueError("Wrong type summary given to MChain!")

        if summary is not None:
            # Check mapping of summary tree and tree set
            if not self.summary.map == self.trees.map:
                try:
                    self.summary.change_mapping(self.trees.map)
                except ValueError as error:
                    raise ValueError(f"{error}\n"
                                     f"The given summary tree and tree set do not fit! "
                                     f"\n(Construction of class MChain failed!)")

    def pwd_matrix(self, csv: bool = False, index="", name=""):
        # todo this should work natively with the name of coupled mchains now, without all these parameters?
        if not os.path.exists(
                f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}.{'csv.gz' if csv else 'npy'}"):
            dm = calc_pw_distances(self.trees)
            if csv:
                np.savetxt(
                    fname=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}.csv.gz",
                    X=dm, delimiter=',', fmt='%i')
            else:
                np.save(
                    file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}.npy",
                    arr=dm)
            return dm
        else:
            if csv:
                return np.genfromtxt(
                    fname=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}.csv.gz",
                    delimiter=',', dtype=int)
            else:
                return np.load(
                    file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}.npy")

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
        upper_i = self.log_data.shape[0] - 1
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
                chain_length = int(self.log_data["Sample"][upper_i]) - int(self.log_data["Sample"][lower_i])
            cur_sampling_interval = int(chain_length / (
                    self.log_data[ess_key][lower_i:(upper_i + 1)].shape[0] - 1 - self.log_data[ess_key][lower_i:(
                    upper_i + 1)].isna().sum()))
            return getattr(ess, f"{ess_method}_ess")(data_list=self.log_data[ess_key][lower_i:(upper_i + 1)].dropna(),
                                                     chain_length=chain_length,
                                                     sampling_interval=cur_sampling_interval)
        else:
            raise ValueError("Not (yet) implemented!")

    def compute_geweke_distances(self, norm: bool = False, add: bool = True, first_range=[0.1, 0.2], last_percent=0.4, index="", name=""):

        new_log_list = gwd.geweke_diagnostic_distances(pw_distances=self.pwd_matrix(index=index, name=name),
                                                       first_range=first_range, last_percent=last_percent)

        if add:
            self.add_new_log_list(new_log_list=new_log_list, col_key=f"Geweke_distances{'_norm' if norm else ''}")
        return new_log_list

    def compute_geweke_focal_tree(self, focal_tree="FM", norm: bool = True, kind="default", add: bool = True,
                                  first_range=[0.1, 0.2], last_percent=0.4):

        new_log_list = gwd.geweke_diagnostic_focal_tree(trees=self.trees, focal_tree=focal_tree, norm=norm, kind=kind,
                                                        first_range=first_range, last_percent=last_percent)

        if add:
            self.add_new_log_list(new_log_list=new_log_list,
                                  col_key=f"Geweke_focal_{focal_tree}_{kind}{'_norm' if norm else ''}")
        return new_log_list

    # todo ideally this should be accessible with the get_ess() funciton and ess_key="pseudo"
    def get_pseudo_ess(self, ess_method="arviz", **kwargs):
        # todo all the testing for upper_i and lower_i but ideally this will be part of the other funciton so this is temporary
        # todo should return the same as RWTY package but also possibly computed with the RNNI distance measure
        if type(ess_method) is str:
            if not hasattr(ess, f"{ess_method}_ess"):
                raise ValueError(f"The given ESS method {ess_method} does not exist!")
        else:
            raise ValueError(ess_method)
        lower_i = 0
        if "lower_i" in kwargs:
            lower_i = kwargs["lower_i"]
        upper_i = len(self.trees)
        if "upper_i" in kwargs:
            upper_i = kwargs["upper_i"]

        if "lower_i" in kwargs or "upper_i" in kwargs:
            if type(lower_i) is not int or type(upper_i) is not int:
                raise ValueError("Wrong type for upper or lower index!")
            if upper_i > self.log_data.shape[0]:
                raise IndexError(f"{upper_i} out of range!")
            if lower_i > upper_i or lower_i < 0 or upper_i < 0:
                raise ValueError("Something went wrong with the given upper and lower index!")

        dist = "rnni"
        if "dist" in kwargs:
            dist = kwargs["dist"]
        sample_range = 10
        if "sample_range" in kwargs:
            sample_range = kwargs["sample_range"]
        chain_length = 1
        if upper_i > 0:
            chain_length = (self.chain_length / (len(self.trees) - 1)) * ((upper_i - 1) - lower_i)
        return ess.pseudo_ess(ess_method=ess_method, tree_set=self.trees[lower_i:upper_i], chain_length=chain_length,
                              sampling_interval=self.chain_length / (len(self.trees) - 1),
                              dist=dist, sample_range=sample_range)

    # todo missing tests
    # the rnni variance log is essentially a running mean plot
    def compute_rnni_variance_log(self, focal_tree_type="tree", norm=False, add=True):
        new_log_list = []
        for i in range(0, len(self.trees)):
            # todo this type checking is not efficient, should be outside of the loop!
            if focal_tree_type == "FM":
                focal_tree = frechet_mean(self.trees[0:i + 1])
            elif focal_tree_type == "tree":
                focal_tree = self.trees[i]
            elif focal_tree_type == "centroid":
                sys.exit("Not yet implemented")
                # focal_tree = frechet_mean(self.trees[0:i + 1])
            else:
                raise ValueError(f"Unknown type {focal_tree_type} given!")
            new_log_list.append(compute_sos_mt(focal_tree, self.trees[0:i + 1], n_cores=None, norm=norm) / (i + 1))

        # adding the computed log value to the log dataframe

        if add:
            self.add_new_log_list(new_log_list=new_log_list, col_key=f"Var{'_norm' if norm else ''}_{focal_tree_type}")
        return new_log_list

    # todo missing test and proper management of the average parameter
    # todo delete the enums things, not really useful in this case
    def compute_new_tree_distance_log(self, average: enums.Average = enums.Average.MEAN, norm: bool = False,
                                      add: bool = True):
        new_log_list = [0]  # initialize as the first iteration is just one tree
        for i in range(1, len(self.trees)):
            new_log_list.append(
                average.function([self.trees[i].fp_distance(self.trees[j], norm=norm) for j in range(0, i)]))
        if add:
            self.add_new_log_list(new_log_list=new_log_list,
                                  col_key=f"Distance{'_norm' if norm else ''}_{average.name}")
        return new_log_list

    def compute_new_tree_summary_distance_log(self, summary="FM", norm=False, add=True):
        new_log_list = [0]  # initialize as the first iteration is just one tree
        for i in range(1, len(self.trees)):
            if summary == "FM":
                new_log_list.append(self.trees[i].fp_distance(frechet_mean(self.trees[0:i]), norm=norm))
            elif summary == "Centroid":
                raise ValueError("Not yet implemented!")
            else:
                raise ValueError(f"Not implemented summary type {summary}!")
        if add:
            self.add_new_log_list(new_log_list=new_log_list, col_key=f"Distance{'_norm' if norm else ''}_{summary}")
        return new_log_list

    # todo missing some proper tests
    def compute_mean_in_chain_deviation(self, norm=False, add=True):
        new_log_list = [1, 1]
        distance_list = {f"{r},{s}": self.trees.fp_distance(r, s, norm=norm) ** 2 for r, s in
                         list(itertools.permutations(range(len(self.trees)), 2))}
        cur_sum = 0
        seen = {}
        for i in range(2, len(self.trees)):
            cur_permutations = list(itertools.permutations(range(i), 2))
            for r, s in cur_permutations:
                if f"{r},{s}" not in seen:
                    cur_sum += distance_list[f"{r},{s}"]
                    seen[f"{r},{s}"] = 1
            new_log_list.append(np.sqrt(cur_sum / len(cur_permutations)))
        if add:
            self.add_new_log_list(new_log_list=new_log_list, col_key=f"In_Chain_deviation{'_norm' if norm else ''}")
        return new_log_list

    # todo missing proper tests
    def compute_ess_traces(self, all=True, partial=[]):
        # computes the ess traces for all the statistic values in the data
        if (not all) and partial:
            # only compute the given list of statistics ess traces
            sys.exit("Currently Feature is not implemented!")
        ess_keys = self.get_key_names()
        method = 'tracerer'  # todo this can be either arviz, coda or tracerer

        for key in ess_keys:
            self.log_data[f"{key}_ess_cum_trace"] = np.nan
            # todo the start of this range/loop could be an argument given
            for i in range(2, self.log_data.shape[0]):
                if not np.isnan(self.log_data[key][i]):
                    self.log_data[f"{key}_ess_cum_trace"][i] = self.get_ess(ess_key=key, ess_method=method, upper_i=i)
        return 0

    # todo missing tests
    # todo this is only applicable if the added list is starting at 0
    #  if it starts with index 5 the sampling interval and everything is messed up
    def add_new_log_list(self, new_log_list, col_key):
        if len(new_log_list) == self.log_data.shape[0]:
            self.log_data[col_key] = new_log_list
        elif len(new_log_list) < self.log_data.shape[0]:
            if self.tree_sampling_interval % self.log_sampling_interval == 0:
                sampling_diff = int(self.tree_sampling_interval / self.log_sampling_interval)
                self.log_data[col_key] = np.nan
                self.log_data[col_key][::sampling_diff] = new_log_list
            else:
                sys.exit("Feature not implemented! Trees logged is not a multiple of the log sampling interval!")
        elif len(new_log_list) > self.log_data.shape[0]:
            sys.exit("Feature not yet implemented! "
                     "The logfile contains less data than the tree file!")

    def norm_column(self, col, normcol, new_key):
        # if value_key not in self.log_data:
        #     raise KeyError(f"Given Value {value_key} does not exist!")
        self.log_data[new_key] = self.log_data[col] / self.log_data[normcol]

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
        _plots._ess_trace_plot(data)
        return 0

    def get_trace_plot(self, value_key='posterior'):
        if value_key not in self.log_data:
            raise KeyError(f"Given Value {value_key} does not exist!")
        _plots._log_trace_plot(self.log_data[value_key][:5])
        return 0

    def get_simmatrix(self, index="", name="", beta=1):
        if not os.path.exists(
                f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}_{'' if beta == 1 else f'{beta}_'}similarity.npy"):

            matrix = self.pwd_matrix()
            similarity = np.exp(-beta * matrix / matrix.std())
            np.save(
                file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_{'' if beta == 1 else f'{beta}_'}similarity.npy",
                arr=similarity)
            return similarity
        else:
            return np.load(
                    file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}_{'' if beta == 1 else f'{beta}_'}similarity.npy")

    def spectral_clustree(self, n_clus=2, beta=1, index="", name="",):
        if not os.path.exists(
                f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_"
                f"{'' if beta == 1 else f'{beta}_'}{n_clus}clustering.npy"):

            clustering = _spectral_clustree(self.get_simmatrix(beta=beta), n_clus=n_clus)
            np.save(
                file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering.npy",
                arr=clustering)
            return clustering
        else:
            return np.load(
                file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering.npy")

    # todo add function that add the ESS traces as columns in the dataframe for the Tracer Visualization
    #  this is easier than writing my own plots

    def write_log_file(self, path):
        # todo check if path is correct
        #  throw warning if file exitsts and stuff like that

        # todo make this a function to calculate all the values and then write a tracer compatible log file
        #  this should also have a parameter accepting columns which will then be included in the log file
        #  default all columns including the original log values from beast logfile

        # todo this should take the working directory into account?

        # This should output a csv file that is compatible with Tracer to visualize all the values
        self.log_data.dropna().to_csv(path, sep="\t", index=False)

    def do_all_the_things(self, out_file, only_norm=True):
        # todo computes all the parameters possible, and then writes a logfile that is compatible with tracer

        # todo should also use parallel processes for each thing!

        # todo write things to the file whenever finished

        # todo should check if the outfile exists and if it does only add the things that are missing!

        # todo Enum for parameter setting instead of the current stuff

        if not only_norm:
            # todo compute all the un normed parameters also
            sys.exit("currently not implemented feature!")

        self.log_data.dropna().to_csv(out_file, sep="\t", index=False)
