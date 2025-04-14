import os
from multiprocessing.managers import Value

import pandas as pd
import numpy as np

from tetres.trees.time_trees import TimeTreeSet
from tetres.judgement._pairwise_distance_matrix import calc_pw_distances
from tetres.clustree.spectral_clustree import _spectral_clustree, spectral_clustree_dm
from tetres.judgement.ess import autocorr_ess, pseudo_ess
from tetres.clustree.bic import bic, plot_bic
from tetres.clustree.silhouette_score import silhouette_score
from tetres.summary.centroid import Centroid
from tetres.judgement.enums import MDS
import warnings

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

        for new_dir in ["data", "plots", "clustering"]:
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

    def pwd_matrix(self, csv: bool = False, name="", rf=False):
        if not os.path.exists(
                f"{self.working_dir}/data/{self.name if name == '' else name}{'_rf' if rf else ''}.{'csv.gz' if csv else 'npy'}"):
            dm = calc_pw_distances(self.trees, rf=rf)
            if csv:
                np.savetxt(
                    fname=f"{self.working_dir}/data/{self.name if name == '' else name}{'_rf' if rf else ''}.csv.gz",
                    X=dm, delimiter=',', fmt='%i')
            else:
                np.save(
                    file=f"{self.working_dir}/data/{self.name if name == '' else name}{'_rf' if rf else ''}.npy",
                    arr=dm)
            return dm
        else:
            if csv:
                return np.genfromtxt(
                    fname=f"{self.working_dir}/data/{self.name if name == '' else name}{'_rf' if rf else ''}.csv.gz",
                    delimiter=',', dtype=int)
            else:
                return np.load(
                    file=f"{self.working_dir}/data/{self.name if name == '' else name}{'_rf' if rf else ''}.npy")

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

    def split_trees_from_clustering(self, k, _overwrite=False):

        # todo change to the get_clustering funciton and also include the type of clustering that is wanted here
        if os.path.exists(f"{self.working_dir}/clustering/sc-{k}-{self.name}.npy"):

            clustering = np.load(f"{self.working_dir}/clustering/sc-{k}-{self.name}.npy")
        else:
            raise ValueError("Clustering has not yet been computed!")

        for k_cluster in range(k):
            cur_treeset = TimeTreeSet()
            cur_treeset.map = self.trees.map
            cur_treeset.trees = [self.trees[index] for index, x in enumerate(clustering) if x == k_cluster]
            if len(cur_treeset) != 0:
                if _overwrite:
                    try:
                        os.remove(f"{self.working_dir}/clustering/trees_k-{k}-c-{k_cluster}-{self.name}.trees")
                    except FileNotFoundError:
                        pass
                if os.path.exists(f"{self.working_dir}/clustering/trees_k-{k}-c-{k_cluster}-{self.name}.trees"):
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

    def get_best_bic_cluster(self, max_cluster=5, local=False):
        # todo add overwrite option and saving this to a file so that it can be read
        #  maybe just save all the bic scores to a file and then create the plot based on that
        best_cluster = 0
        best_bic = None
        prev_mnd = None
        # todo soon to be arguments ?
        beta = 1

        local_norm = False
        if local:
            local_norm = np.max(self.pwd_matrix())

        for cur_cluster in range(max_cluster):

            summaries = self.get_cluster_centroids(k=cur_cluster+1, beta=beta)
            clustering = self.get_clustering(k=cur_cluster+1)

            cur_bic, mnd_cluster = bic(treeset=self.trees,
                                       clustering=clustering,
                                       summaries=summaries,
                                       local_norm=local_norm)

            if prev_mnd is not None and sum(mnd_cluster.values()) > prev_mnd:
                warnings.warn(f"MND Constraint reached at {cur_cluster+1} clustering")
                return best_cluster
            prev_mnd = sum(mnd_cluster.values())
            if best_bic is None or cur_bic < best_bic:
                best_cluster = cur_cluster + 1
                best_bic = cur_bic
        return best_cluster

    def get_best_silhouette_cluster(self, max_cluster=5, local_norm=False):
        best_cluster = 0
        best_sil = None
        # todo rework with new funcitons
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

    def get_clustering(self, k, _overwrite=False, beta=1):

        cluster_type = f"sc{'' if beta == 1 else f'-b{beta}'}"
        # todo should be a list of cluster_type [sc, ...] Future feature
        # todo beta should be part of kwargs as only needed for similarity matrix calculation for spectral clustering

        if k == 1:
            return np.zeros(len(self.trees), dtype=int)
        elif k > 1:
            if _overwrite:
                try:
                    os.remove(f"{os.path.join(self.working_dir, 'clustering')}/{cluster_type}-{k}-{self.name}.npy")
                except FileNotFoundError:
                    pass

            if os.path.exists(f"{os.path.join(self.working_dir, 'clustering')}/{cluster_type}-{k}-{self.name}.npy"):
                clustering = np.load(os.path.join(self.working_dir, 'clustering', f"{cluster_type}-{k}-{self.name}.npy"))
            else:
                if cluster_type == f"sc{'' if beta == 1 else f'-b{beta}'}":
                    clustering = _spectral_clustree(self.similarity_matrix(beta=beta), n_clus=k)
                    np.save(file=f"{os.path.join(self.working_dir, 'clustering')}/{cluster_type}-{k}-{self.name}", arr=clustering)
            return clustering
        else:
            raise ValueError(f"Value {k} not supported!")

    def get_cluster_centroids(self, k, beta=1, _overwrite=False, _best_of_range = 10):
        # todo beta as kwargs parameter in the future
        cluster_type = f"sc{'' if beta == 1 else f'-b{beta}'}"

        output = []
        for cur_cluster in range(k):
            cur_centroid_path = os.path.join(self.working_dir, "clustering", f"{cluster_type}-{k}{f'-k{cur_cluster+1}' if k != 1 else ''}-{self.name}.tree")
            if _overwrite:
                try:
                    os.remove(cur_centroid_path)
                except FileNotFoundError:
                    pass

            # Creating a TimeTreeSet for the k-th cluster
            cur_treeset = TimeTreeSet()
            cur_treeset.map = self.trees.map
            cur_treeset.trees = [self.trees[index] for index, x in
                                         enumerate(self.get_clustering(k=cur_cluster + 1)) if x == cur_cluster]
            if len(cur_treeset) == 0:
                output.append(None)
            else:
                if os.path.exists(cur_centroid_path):
                    # Assumes that the map is always the same
                    output.append(TimeTreeSet(file=cur_centroid_path)[0])
                else:
                    centroid = Centroid()
                    cen, sos = centroid.compute_centroid(cur_treeset)
                    for _ in range(_best_of_range):
                        new_cen, new_sos = centroid.compute_centroid(cur_treeset)
                        if new_sos < sos:
                            cen, sos = new_cen, new_sos
                    cen.write_nexus(cur_treeset.map,
                                    cur_centroid_path,
                                    name=f"Cen-{k}-k{cur_cluster}-{self.name}")
                    output.append(cen)
        return output

    def get_mds_coords(self, mds_type = 'tsne', dim=2, dist='rf'):
        # todo will need to add more possible parameters, kwargs...
        # todo sort out how to properly hande dist= make it an enum and generally expandable for different distances...

        # todo this requires a proper test
        if not isinstance(mds_type, MDS):
            if mds_type in MDS:
                mds_type = MDS[mds_type]
            else:
                raise TypeError('mds_type invalid.')

        mds_coords_filename = 'BLBALBALBA'  # f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity.npy"

        if not os.path.exists(mds_coords_filename):
            # need to first calculate the coordinates with the correct funciton
            pwd_matrix = self.pwd_matrix(rf=False)  # todo make a dist enum and change the respective function for this
            if mds_type == MDS.TSNE:
                coords, kl_divergence = mds_type.value(pwd_matrix=pwd_matrix, dim=dim)
            elif mds_type == MDS.PYTHON:
                raise NotImplementedError('Python MDS not implemented yet.')
            elif mds_type == MDS.R:
                raise NotImplementedError('R isoMDS not implemented yet.')  # todo there is also cmdscale?
            # todo could add spectral embedding, very closely related to the spectral clustering
            else:
                raise ValueError(f'Unsupported mds_type argument provided {mds_type}')
            np.save(file=mds_coords_filename, arr=coords)
        else:
            # have to read the already computed coords
            coords = np.load(file=mds_coords_filename)
        return coords


        # todo support for different types of mds coords, saving them and reading them
        #  different distances
        #  dimensionality is also an option

        # # computing or loading the tsne coordinates
        # # todo maybe change this to do a recomputation if wanted and overwriting of the old if better fit
        #
        # # todo this should just be a funciton that computes tsne, all the questions should be part of the mchain function
        # coords_output = f"{mchain.working_dir}/data/{mchain.name}_tsne{'_rf' if rf else ''}_{dim}d_coords.npy"
        # if not os.path.exists(coords_output):
        #     pwd_matrix = mchain.pwd_matrix(rf=rf)
        #     coords, kl_divergence = _tsne_coords_from_pwd(pwd_matrix=pwd_matrix, dim=dim)
        #     # todo do something with the kl_divergence value
        #     np.save(file=coords_output, arr=coords)
        #     return coords
        # else:
        #     return np.load(coords_output)

    def similarity_matrix(self, beta=1):
        if not os.path.exists(
                f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity.npy"):
            matrix = self.pwd_matrix()
            if np.allclose(matrix, np.triu(matrix)):
                # matrix is upper triangular, has to be changed to a full matrix for this implementation
                matrix += np.transpose(matrix)
            similarity = np.exp(-beta * matrix / matrix.std())
            np.save(
                file=f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity.npy",
                arr=similarity)
            return similarity
        else:
            return np.load(
                file=f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity.npy")

    def evaluate_clustering(self, kind="bic", _overwrite_clustering=False, add_random=True, local=False):

        # todo parameters for overwrite currently not used

        # todo kind should be an enum
        if kind not in ["silhouette", "bic"]:
            raise ValueError(f"Kind {kind} not supported, choose either 'Silhouette' or 'BIC'.")

        # todo if _overwrite_plot delete plot and silhouette score
        # todo if _recalculate_clsutering delete clustering and also the plot and save file for the scores

        if kind == "bic":
            # todo kwargs with add_random and local norm, max_cluster
            max_cluster = 5
            beta = 1
            if local:
                # Set local norm to be the max distance observed instead of diameter of the space
                local_norm = np.max(self.pwd_matrix())
            else:
                local_norm = False

            clusterings = {k+1: self.get_clustering(k+1, beta=beta, _overwrite=_overwrite_clustering) for k in range(max_cluster)}
            summaries_dict = {k+1: self.get_cluster_centroids(k=k+1, _overwrite=_overwrite_clustering) for k in range(max_cluster)}

            # _overwrite=True will recalculate the clustering too, as this is being passed down to the bic() function
            plot_bic(treeset=self.trees,
                     clusterings=clusterings,
                     summaries_dict=summaries_dict,
                     max_cluster=max_cluster,
                     local_norm=local_norm,
                     add_random=add_random,
                     file_name=os.path.join(self.working_dir, "plots", f"BIC_{'random_' if add_random else ''}{'local_' if local_norm else ''}{self.name}.pdf"))

