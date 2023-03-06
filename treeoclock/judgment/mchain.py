import sys
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import itertools
import linecache
import re
import warnings

from treeoclock.trees.time_trees import TimeTreeSet
from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary.frechet_mean import frechet_mean
from treeoclock.trees._converter import ctree_to_ete3
from treeoclock.judgment._pairwise_distance_matrix import calc_pw_distances, calc_pw_distances_two_sets
from treeoclock.judgment import _gelman_rubin_diag as grd
from treeoclock.judgment import _cladesetcomparator as csc
from treeoclock.clustree.spectral_clustree import _spectral_clustree
from treeoclock.summary.centroid import Centroid
from treeoclock.summary.annotate_centroid import annotate_centroid
from treeoclock.judgment.ess import autocorr_ess, pseudo_ess
from treeoclock.judgment._discrete_cladesetcomparator import discrete_cladeset_comparator

from rpy2.robjects.packages import importr
phytools = importr("phytools")
phangorn = importr("phangorn")


class coupled_MChains():
    def __init__(self, m_MChains, trees, log_files, working_dir, name="cMC"):
        # todo currently still missing the summary parameter of MChain, but its not really used at this point
        self.m_MChains = m_MChains
        self.name = name

        # if type(trees) != type(log_files):
        #     raise ValueError("trees and logfiles should be given with the same data type!")

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
            self.tree_files = trees.copy()
            self.log_files = log_files.copy()
            for i in range(self.m_MChains):
                self.MChain_list.append(MChain(trees=trees[i], log_file=log_files[i],
                                               working_dir=working_dir, name=f"{self.name}_{i}"))

        # todo check that the mappings are all the same for every chain if not change them to one!

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

    def pwd_matrix(self, index1, index2=None, csv: bool = False, rf: bool = False):
        if type(index1) is not int:
            raise ValueError("Unrecognized index type!")
        if index1 >= self.m_MChains:
            raise IndexError("Given Index out of range!")

        if index2 is None:
            return self.MChain_list[index1].pwd_matrix(index=index1, name=self.name, csv=csv, rf=rf)
        else:
            if type(index2) is not int:
                raise ValueError("Unrecognized index type!")
            if index1 == index2:
                raise ValueError("Indeces are equal, only call with one index for pairwise distances of one set!")
            if index2 >= self.m_MChains:
                raise IndexError("Given Index out of range!")
            if index2 < index1:
                index1, index2 = index2, index1
            if not os.path.exists(f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.{'csv.gz' if csv else 'npy'}"):
                dm = calc_pw_distances_two_sets(self.MChain_list[index1].trees, self.MChain_list[index2].trees, rf=rf)
                if csv:
                    np.savetxt(fname=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.csv.gz", X=dm, delimiter=',',
                               fmt='%i')
                else:
                    np.save(file=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.npy", arr=dm)
                return dm
            else:
                if csv:
                    return np.genfromtxt(fname=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.csv.gz",
                                         delimiter=',', dtype=int)
                else:
                    return np.load(f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.npy")

    # todo testing required
    def pwd_matrix_all(self, rf: bool = False):
        if not os.path.exists(f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy"):
            combined_matrix = np.array([])
            for i in range(self.m_MChains):
                cur_row = np.array([])
                for j in range(self.m_MChains):
                    # print(i, j)
                    if i < j:
                        cur_row = np.concatenate((cur_row, self.pwd_matrix(i, j, rf=rf)),
                                                 axis=1) if cur_row.size else self.pwd_matrix(i, j, rf=rf)
                    elif i > j:
                        cur_row = np.concatenate((cur_row, np.zeros((len(self[i].trees), len(self[j].trees)))),
                                                 axis=1) if cur_row.size else np.zeros((len(self[i].trees), len(self[j].trees)))
                    elif i == j:
                        cur_row = np.concatenate((cur_row, self.pwd_matrix(i, rf=rf)),
                                                 axis=1) if cur_row.size else self.pwd_matrix(i, rf=rf)
                    # print(cur_row.shape)
                combined_matrix = np.concatenate((combined_matrix, cur_row), axis=0) if combined_matrix.size else cur_row
            np.save(file=f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy", arr=combined_matrix)
            return combined_matrix
        else:
            return np.load(f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy")

    def _extract_cutoff(self, i, start, end, ess, subsample, compare_to, _overwrite=False):
        # todo should be its own script, also rename all the parameters given
        try:
            os.mkdir(f"{self.working_dir}/cutoff_files")
        except FileExistsError:
            pass

        tree_file = f"{self.working_dir}/cutoff_files/{self[i].name}_{self[compare_to].name}{'' if ess == 0 else f'_ess-{ess}'}{f'_subsample-{subsample}' if subsample else ''}.trees"
        log_file = f"{self.working_dir}/cutoff_files/{self[i].name}_{self[compare_to].name}{'' if ess == 0 else f'_ess-{ess}'}{f'_subsample-{subsample}' if subsample else ''}.log"

        if _overwrite:
            try:
                os.remove(tree_file)
            except FileNotFoundError:
                pass
            try:
                os.remove(log_file)
            except FileNotFoundError:
                pass

        if self.log_files[i] is not None:
            try:
                with open(log_file, "x") as f:
                    f.write("\t".join(v for v in list(self[i].log_data.keys())))
                    f.write("\n")
            except FileExistsError:
                return 0  # assuming this was already done hence skipping this function

        try:
            with open(tree_file, "x") as f:
                f.write(f"#NEXUS\n\nBegin taxa;\n\tDimensions ntax={self[i].trees[i].ctree.num_leaves};\n\t\tTaxlabels\n")
                for taxa in range(1, self[i].trees[0].ctree.num_leaves + 1):
                    f.write(f"\t\t\t{self[i].trees.map[taxa]}\n")
                f.write("\t\t\t;\nEnd;\nBegin trees;\n\tTranslate\n")
                for taxa in range(1, self[i].trees[0].ctree.num_leaves):
                    f.write(f"\t\t\t{taxa} {self[i].trees.map[taxa]},\n")
                f.write(
                    f"\t\t\t{self[i].trees[0].ctree.num_leaves} {self[i].trees.map[self[i].trees[0].ctree.num_leaves]}\n")
                f.write(";\n")
        except FileExistsError:
            return 0  # assuming it was already done, therefore end function

        sample = 1
        chain_treefile = f"{self.working_dir}/{self.tree_files[i]}"
        for index in range(start, end+1):
            re_tree = re.compile("\t?tree .*=? (.*$)", flags=re.I | re.MULTILINE)
            offset = 10 + (2 * self[i].trees[0].ctree.num_leaves)

            line = index + 1 + offset
            cur_tree = linecache.getline(chain_treefile, line)
            cur_tree = f'{re.split(re_tree, cur_tree)[1][:re.split(re_tree, cur_tree)[1].rfind(")") + 1]};'
            with open(tree_file, "a") as tf:
                tf.write(f"tree STATE_{sample} = {cur_tree}\n")

            if self.log_files[i] is not None:
                cur_values = [v for v in self[i].log_data.iloc[index]]
                cur_values[0] = sample
                with open(log_file, "a") as lf:
                    lf.write("\t".join([str(v) for v in cur_values]))
                    lf.write("\n")
            sample += 1
        with open(tree_file, "a") as tf:
            tf.write("End;")

    def gelman_rubin_all_chains_density_plot(self, samples: int = 100):
        if len(self) < 2:
            raise ValueError("Gelman Rubin only possible with multiple Chains!")
        return grd.gelman_rubin_all_chains_density_plot(self, samples=samples)

    def gelman_rubin_parameter_choice_plot(self, i, j, _subsampling=False):
        if len(self) < 2:
            raise TypeError("Gelman Rubin only possible with multiple Chains!")
        return grd.gelman_rubin_parameter_choice_plot(self, i, j, _subsampling=_subsampling)

    def gelman_rubin_cut(self, i, j, smoothing=0.5, ess_threshold=0, pseudo_ess_range=100, _overwrite=False, smoothing_average="mean", _subsampling=False):
        # First if overwrite delete previous computation file
        if _overwrite:
            try:
                os.remove(
                    f"{self.working_dir}/data/{self.name}_{i}_{j}_gelman_rubin_cutoff{f'_subsampling-{_subsampling}' if _subsampling else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}_smoothing-{smoothing}_{smoothing_average}")
            except FileNotFoundError:
                pass
        # if the file still exists here then we do not need to compute anything
        if os.path.exists(f"{self.working_dir}/data/{self.name}_{i}_{j}_gelman_rubin_cutoff{f'_subsampling-{_subsampling}' if _subsampling else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}_smoothing-{smoothing}_{smoothing_average}"):
            with open(
                    f"{self.working_dir}/data/{self.name}_{i}_{j}_gelman_rubin_cutoff{f'_subsampling-{_subsampling}' if _subsampling else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}_smoothing-{smoothing}_{smoothing_average}",
                    "r") as file:
                cut_start = int(file.readline())
                cut_end = int(file.readline())
            return cut_start, cut_end

        cut_start, cut_end = grd.gelman_rubin_cut(self, i=i, j=j,
                                                  smoothing=smoothing,
                                                  ess_threshold=ess_threshold,
                                                  pseudo_ess_range=pseudo_ess_range,
                                                  smoothing_average=smoothing_average,
                                                  _subsampling=_subsampling)
        # Write the cutoff boundaries to a file, if it already exists skip this part
        if _subsampling:
            for _ in range(_subsampling):
                cut_start, cut_end = cut_start * 2, cut_end * 2
        try:
            with open(
                    f"{self.working_dir}/data/{self.name}_{i}_{j}_gelman_rubin_cutoff{f'_subsampling-{_subsampling}' if _subsampling else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}_smoothing-{smoothing}_{smoothing_average}",
                    "x") as f:
                f.write(f"{cut_start}\n{cut_end}")
        except FileExistsError:
            raise Exception("Something went wrong!")
        return cut_start, cut_end

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
        # todo add the compare chronogram on the lower diag of the plot?

    def cen_for_each_chain(self, repeat=5):
        for _ in range(repeat):
            for chain in self.MChain_list:
                cen = Centroid(variation="inc_sub")
                cen_tree, sos = cen.compute_centroid(chain.trees)
                cen_tree = annotate_centroid(cen_tree, chain.trees)
                try:
                    file = open(f"{chain.working_dir}/data/{chain.name}_cen_sos.log", "r")
                    cur_sos = int(file.readline())
                    file.close()
                except FileNotFoundError:
                    # Setting cur_sos so that the folowing if will be Treu
                    cur_sos = sos +1
                if sos < cur_sos:
                    print("found better centroid!")
                    # Removing the old sos file if exists
                    try:
                        os.remove(f"{chain.working_dir}/data/{chain.name}_cen_sos.log")
                    except FileNotFoundError:
                        pass
                    # Writing the better sos value to the file
                    file = open(f"{chain.working_dir}/data/{chain.name}_cen_sos.log", "x")
                    file.write(f"{sos}")
                    file.close()
                    # deleting the old centroid tree file
                    try:
                        os.remove(f"{chain.working_dir}/{chain.name}_cen.tree")
                    except FileNotFoundError:
                        pass
                    # writing the better centroid tree to the file
                    cen_tree.write_nexus(chain.trees.map, file_name=f"{chain.working_dir}/{chain.name}_cen.tree", name=f"Cen_{chain.name}")

    def clade_set_comparison(self, i, j, plot=True, burnin=0):
        return discrete_cladeset_comparator(tree_set_i=self[i].trees, tree_set_j=self[j].trees, plot=plot, burnin=burnin,
                                            file=f"{self.working_dir}/plots/{self.name}_discrete_cc_{i}_{j}_burn-{burnin}.png")


class MChain:
    def __init__(self, working_dir, trees, log_file=None, summary=None, name: str = "MC"):
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
            self.trees = TimeTreeSet(f"{self.working_dir}/{trees}")  # todo use os.path.join
        else:
            raise ValueError(f"{type(trees)} of trees not recognized!")

        # Reading the log_file
        if log_file is None:
            # if no logfile set chain length and tree sampling interval
            self.chain_length = len(self.trees) - 1
            self.tree_sampling_interval = 1
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
            elif log_file is not None:
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

    def __len__(self):
        # returns number of tree samples at this point!
        return len(self.trees)

    def pwd_matrix(self, csv: bool = False, index="", name="", rf=False):
        # todo this should work natively with the name of coupled mchains now, without all these parameters?
        if not os.path.exists(
                f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}{'_rf' if rf else ''}.{'csv.gz' if csv else 'npy'}"):
            dm = calc_pw_distances(self.trees, rf=rf)
            if csv:
                np.savetxt(
                    fname=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}{'_rf' if rf else ''}.csv.gz",
                    X=dm, delimiter=',', fmt='%i')
            else:
                np.save(
                    file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}{'_rf' if rf else ''}.npy",
                    arr=dm)
            return dm
        else:
            if csv:
                return np.genfromtxt(
                    fname=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}{'_rf' if rf else ''}.csv.gz",
                    delimiter=',', dtype=int)
            else:
                return np.load(
                    file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index!='' else ''}{'_rf' if rf else ''}.npy")

    def get_key_names(self):
        return list(self.log_data.columns)[1:]

    def get_ess(self, ess_key="posterior", **kwargs):
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
            return autocorr_ess(data_list=list(self.log_data[ess_key][lower_i:(upper_i + 1)].dropna()))
        else:
            raise ValueError("Not (yet) implemented!")

    def split_trees_from_clustering(self, k):
        if os.path.exists(f"{self.working_dir}/clustering/sc-{k}-{self.name}.npy"):
            clustering = np.load(f"{self.working_dir}/clustering/sc-{k}-{self.name}.npy")
        else:
            raise ValueError("Clustering has not yet been computed!")

        for k_cluster in range(k):
            cur_treeset = TimeTreeSet()
            cur_treeset.map = self.trees.map
            cur_treeset.trees = [self.trees[index] for index, x in enumerate(clustering) if x == k_cluster]
            if len(cur_treeset) != 0:
                if os.path.exists(f"{self.working_dir}/clustering/trees_k-{k}-c-{k_cluster}.trees"):
                    raise ValueError("Tree file already exists!")
                cur_treeset.write_nexus(file_name=f"{self.working_dir}/clustering/trees_k-{k}-c-{k_cluster}.trees")

    def get_pseudo_ess(self, **kwargs):
        lower_i = 0
        if "lower_i" in kwargs:
            lower_i = kwargs["lower_i"]
        upper_i = len(self.trees)
        if "upper_i" in kwargs:
            upper_i = kwargs["upper_i"]

        if "lower_i" in kwargs or "upper_i" in kwargs:
            if type(lower_i) is not int or type(upper_i) is not int:
                raise ValueError("Wrong type for upper or lower index!")
            if upper_i > len(self.trees):
                raise IndexError(f"{upper_i} out of range!")
            if lower_i > upper_i or lower_i < 0 or upper_i < 0:
                raise ValueError("Something went wrong with the given upper and lower index!")

        dist = "rnni"
        if "dist" in kwargs:
            dist = kwargs["dist"]
        sample_range = 10
        if "sample_range" in kwargs:
            sample_range = kwargs["sample_range"]
        no_zero = False
        if "no_zero" in kwargs:
            no_zero = kwargs["no_zero"]
        # chain_lngth = 1
        # if upper_i > 0:
        #     chain_leength = (self.chain_length / (len(self.trees) - 1)) * ((upper_i - 1) - lower_i)
        return pseudo_ess(tree_set=self.trees[lower_i:upper_i],
                              dist=dist, sample_range=sample_range, no_zero=no_zero)

    def simmilarity_matrix(self, index="", name="", beta=1):

        # todo adapt to new folder clustering/ and integrate what I did for BIC here!

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

        # todo adapt to new folder clustering/ and integrate what I did for BIC here!

        if not os.path.exists(
                f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_"
                f"{'' if beta == 1 else f'{beta}_'}{n_clus}clustering.npy"):

            clustering = _spectral_clustree(self.simmilarity_matrix(beta=beta), n_clus=n_clus)
            np.save(
                file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering.npy",
                arr=clustering)
            return clustering
        else:
            return np.load(
                file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering.npy")

    def evaluate_clusterin(self, kind="Silhouette"):
        if kind not in ["Silhouette", "BIC"]:
            raise ValueError(f"Kind {kind} not supported, choose either 'Silhouette' or 'BIC'.")
        # todo implementation missing
        return 0
