import os
import pandas as pd
import numpy as np

from tetres.trees.time_trees import TimeTreeSet
from tetres.judgment._pairwise_distance_matrix import calc_pw_distances
from tetres.clustree.spectral_clustree import _spectral_clustree
from tetres.judgment.ess import autocorr_ess, pseudo_ess
from tetres.clustree.bic import BIC
from tetres.clustree.silhouette_score import silhouette_score


class Chain:
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

        for new_dir in ["data", "plots"]:
            try:
                os.mkdir(f"{self.working_dir}/{new_dir}")
            except FileExistsError:
                pass

        # Setting trees
        if type(trees) is TimeTreeSet:
            self.trees = trees
        elif type(trees) is str:
            self.trees = TimeTreeSet(os.path.join(self.working_dir, trees))
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
                cur_treeset.write_nexus(file_name=f"{self.working_dir}/clustering/trees_k-{k}-c-{k_cluster}-{self.name}.trees")

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
        return pseudo_ess(tree_set=self.trees[lower_i:upper_i],
                          dist=dist, sample_range=sample_range, no_zero=no_zero)

    def simmilarity_matrix(self, index="", name="", beta=1):

        # todo adapt to new folder clustering/ and integrate what I did for BIC here!

        # todo this matrix should also be in the clustering folder

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
        raise NotImplemented("Currently WIP")

    def evaluate_clustering(self, kind="silhouette", plot=True, _overwrite_plot=False, _overwrite_clustering=False):
        if kind not in ["silhouette", "bic"]:
            raise ValueError(f"Kind {kind} not supported, choose either 'Silhouette' or 'BIC'.")

        # todo if _overwrite_plot delete plot and silhouette score
        # todo if _recalculate_clsutering delete clustering and also the plot and save file for the scores


        # todo change the BIC and silhouette funcitons to take a clustering, reading should take place here instead!

        # todo file identifyer is kind_....
        # todo this function should return the suggeted number of clusters for a given chain, i.e. save the scores to a file and read from that file

        # todo implementation missing
        raise NotImplemented("Currently WIP")

    def get_best_bic_cluster(self, max_cluster=5, local_norm=False):
        # todo add overwrite option and saving this to a file so that it can be read
        #  maybe just save all the bic scores to a file and then create the plot based on that
        best_cluster = 0
        best_bic = None
        for cur_cluster in range(max_cluster):
            bic, mnd_cluster = BIC(treeset=self.trees, matrix=self.pwd_matrix(),
                                k=cur_cluster+1, local_norm=local_norm,
                                working_folder=self.working_dir,
                                _overwrite=False,
                                random_shuffle=False,
                                chain_id=self.name
                                )
            # todo at some point add check of mnd where the total sum should decrease and
            #  an increase indicates that this is not useful any longer
            if best_bic is None or bic < best_bic:
                    best_cluster = cur_cluster + 1
                    best_bic = bic
        return best_cluster

    def get_best_silhouette_cluster(self, max_cluster=5, local_norm=False):
        best_cluster = 0
        best_sil = None
        for cur_cluster in range(2, max_cluster+1):
            cur_s = silhouette_score(matrix=self.pwd_matrix(),
                                     k=cur_cluster,
                                     local_norm=local_norm,
                                     working_folder=self.working_dir,
                                     random_shuffle=False,
                                     chain_id=self.name)

            if best_sil is None or cur_s < best_sil:
                best_cluster = cur_cluster + 1
                best_sil = cur_s
        return best_cluster

    def spectral_clustree(self, n_clus=2, beta=1):
        raise NotImplementedError("Currently not supported and a WIP")
        # todo adapt to new folder clustering/ and integrate what I did for BIC here!
        try:
            os.mkdir(os.path.join(self.working_dir, "clustering"))
        except FileExistsError:
            raise ValueError("Didn't work?!")

        if not os.path.exists(
                f"{self.working_dir}/clustering/sc-{n_clus}-{self.name}{'' if beta == 1 else f'-{beta}'}.npy"):

            clustering = _spectral_clustree(self.simmilarity_matrix(beta=beta), n_clus=n_clus)
            np.save(
                file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering.npy",
                arr=clustering)
            return clustering
        else:
            return np.load(
                file=f"{self.working_dir}/data/{self.name if name == '' else name}{f'_{index}' if index != '' else ''}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering.npy")