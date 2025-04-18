import os
import warnings

import pandas as pd
import numpy as np

from tetres.trees.time_trees import TimeTreeSet
from tetres.judgement._pairwise_distance_matrix import calc_pw_distances
from tetres.utils.decorators import validate_literal_args
from tetres.visualize.mds_coord_compuation import _tsne_coords_from_pwd
from tetres.clustree.spectral_clustree import _spectral_clustree
from tetres.judgement.ess import autocorr_ess, pseudo_ess
from tetres.clustree.bic import bic, plot_bic
from tetres.clustree.silhouette_score import silhouette_score
from tetres.summary.centroid import Centroid
from tetres.utils.literals import _DIST, _MDS_TYPES
from tetres.visualize.plot_config import PlotOptions
from tetres.visualize.plot_coords import plot_coords


# todo this should be moved to the visualize module...



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
            self.log_data = pd.read_csv(f"{self.working_dir}/{log_file}", header=0, sep=r"\s+",
                                        comment="#")
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
            return autocorr_ess(
                data_list=list(self.log_data[ess_key][lower_i:(upper_i + 1)].dropna()))
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
            cur_treeset.trees = [self.trees[index] for index, x in enumerate(clustering) if
                                 x == k_cluster]
            if len(cur_treeset) != 0:
                if _overwrite:
                    try:
                        os.remove(
                            f"{self.working_dir}/clustering/trees_k-{k}-c-{k_cluster}-{self.name}.trees")
                    except FileNotFoundError:
                        pass
                if os.path.exists(
                        f"{self.working_dir}/clustering/trees_k-{k}-c-{k_cluster}-{self.name}.trees"):
                    raise ValueError("Tree file already exists!")
                cur_treeset.write_nexus(
                    file_name=f"{self.working_dir}/clustering/trees_k-{k}-c-{k_cluster}-{self.name}.trees")

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
        if "dist_type" in kwargs:
            dist = kwargs["dist_type"]
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

            summaries = self.get_cluster_centroids(k=cur_cluster + 1, beta=beta)
            clustering = self.get_clustering(k=cur_cluster + 1)

            cur_bic, mnd_cluster = bic(treeset=self.trees,
                                       clustering=clustering,
                                       summaries=summaries,
                                       local_norm=local_norm)

            if prev_mnd is not None and sum(mnd_cluster.values()) > prev_mnd:
                warnings.warn(f"MND Constraint reached at {cur_cluster + 1} clustering")
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
        for cur_cluster in range(2, max_cluster + 1):
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
                    os.remove(
                        f"{os.path.join(self.working_dir, 'clustering')}/{cluster_type}-{k}-{self.name}.npy")
                except FileNotFoundError:
                    pass

            if os.path.exists(
                    f"{os.path.join(self.working_dir, 'clustering')}/{cluster_type}-{k}-{self.name}.npy"):
                clustering = np.load(os.path.join(self.working_dir, 'clustering',
                                                  f"{cluster_type}-{k}-{self.name}.npy"))
            else:
                if cluster_type == f"sc{'' if beta == 1 else f'-b{beta}'}":
                    clustering = _spectral_clustree(self.similarity_matrix(beta=beta), n_clus=k)
                    np.save(
                        file=f"{os.path.join(self.working_dir, 'clustering')}/{cluster_type}-{k}-{self.name}",
                        arr=clustering)
            return clustering
        else:
            raise ValueError(f"Value {k} not supported!")

    def get_cluster_centroids(self, k, beta=1, _overwrite=False, _best_of_range=10):
        # todo beta as kwargs parameter in the future
        cluster_type = f"sc{'' if beta == 1 else f'-b{beta}'}"

        output = []
        for cur_cluster in range(k):
            cur_centroid_path = os.path.join(self.working_dir, "clustering",
                                             f"{cluster_type}-{k}{f'-k{cur_cluster + 1}' if k != 1 else ''}-{self.name}.tree")
            if _overwrite:
                try:
                    os.remove(cur_centroid_path)
                except FileNotFoundError:
                    pass

            # Creating a TimeTreeSet for the k-th cluster
            cur_treeset = TimeTreeSet()
            cur_treeset.map = self.trees.map
            cur_treeset.trees = [self.trees[index] for index, x in
                                 enumerate(self.get_clustering(k=cur_cluster + 1)) if
                                 x == cur_cluster]
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

    @validate_literal_args(mds_type=_MDS_TYPES, dist_type=_DIST)
    def get_mds_coords(self, mds_type: _MDS_TYPES = 'tsne',
                       dim: int = 2, dist_type: _DIST = 'rnni',
                       _overwrite: bool = False) -> np.ndarray:

        warnings.warn("Currently not fully implemented.. Will only compute RNNI TSNE 2dim MDS.",
                      stacklevel=3)

        # todo will need to add more possible parameters, kwargs... beta for spectral

        mds_coords_filename = os.path.join(self.working_dir, "data",
                                           f"{self.name}_{mds_type}-{dist_type}.npy")

        if not os.path.exists(mds_coords_filename) or _overwrite:
            # need to first calculate the coordinates with the correct function

            # todo no really using the mds_type atm, need to implement that
            coords = _tsne_coords_from_pwd(self.pwd_matrix(rf=False), dim=dim)
            np.save(file=mds_coords_filename, arr=coords)
        else:
            # have to read the already computed coords
            coords = np.load(file=mds_coords_filename)
        return coords

    def similarity_matrix(self, beta=1):
        if not os.path.exists(
                f"{self.working_dir}/data/{self.name}_"
                f"{'' if beta == 1 else f'{beta}_'}similarity.npy"):
            matrix = self.pwd_matrix()
            if np.allclose(matrix, np.triu(matrix)):
                # matrix is upper triangular,
                # has to be changed to a full matrix for this implementation
                matrix += np.transpose(matrix)
            similarity = np.exp(-beta * matrix / matrix.std())
            np.save(
                file=f"{self.working_dir}/data/{self.name}_"
                     f"{'' if beta == 1 else f'{beta}_'}similarity.npy",
                arr=similarity)
            return similarity
        else:
            return np.load(
                file=f"{self.working_dir}/data/{self.name}_"
                     f"{'' if beta == 1 else f'{beta}_'}similarity.npy")

    def evaluate_clustering(self, kind="bic", _overwrite_clustering=False, add_random=True,
                            local=False):

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

            clusterings = {
                k + 1: self.get_clustering(k + 1, beta=beta, _overwrite=_overwrite_clustering) for k
                in range(max_cluster)}
            summaries_dict = {
                k + 1: self.get_cluster_centroids(k=k + 1, _overwrite=_overwrite_clustering) for k
                in range(max_cluster)}

            # _overwrite=True will recalculate the clustering too, as this is being passed down to the bic() function
            plot_bic(treeset=self.trees,
                     clusterings=clusterings,
                     summaries_dict=summaries_dict,
                     max_cluster=max_cluster,
                     local_norm=local_norm,
                     add_random=add_random,
                     file_name=os.path.join(self.working_dir, "plots",
                                            f"BIC_{'random_' if add_random else ''}{'local_' if local_norm else ''}{self.name}.pdf"))

    def plot_mds(self, mds_type: _MDS_TYPES = 'tsne', dim: int = 2,
                 dist_type: _DIST = 'rnni') -> None:
        """
        (WIP) Plotting simple MDS of a chain without any clustering.

        :param mds_type: Type of MDS to plot
        :param dim: Dimension to plot to (only dim=2 supported atm)
        :param dist_type: Type of tree distance to use (RNNI or RF)
        :return: None
        """

        # todo the function should take as input the PlotOptions options argument
        #  then here only set filename if not already given
        if dim != 2:
            raise ValueError("Unsupported dimension for MDS -- only supports dim=2 at the moment.")

        mds_plot_file = os.path.join(self.working_dir, "plots", f"{self.name}_{mds_type}-"
                                                                f"{dist_type}.pdf")

        coords = self.get_mds_coords(mds_type=mds_type, dim=dim, dist_type=dist_type)

        # todo title and more information to be put into the plot.
        plot_coords(coords=coords, dim=dim, options=PlotOptions(filename=mds_plot_file))
