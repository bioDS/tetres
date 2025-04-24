import os
import warnings

import pandas as pd
import numpy as np

from tetres.trees.time_trees import TimeTreeSet
from tetres.judgement._pairwise_distance_matrix import calc_pw_distances
from tetres.utils.decorators import validate_literal_args
from tetres.visualize.mds_coord_compuation import _tsne_coords_from_pwd
from tetres.clustree.spectral_clustree import spectral_clustree_dm
from tetres.judgement.ess import calc_ess, pseudo_ess
from tetres.clustree.bic import bic, plot_bic
from tetres.clustree.silhouette_score import silhouette_score
from tetres.utils.literals import DIST, MDS_TYPES, CLUSTERING_TYPE
from tetres.visualize.plot_config import PlotOptions
from tetres.visualize.plot_coords import plot_coords


class Chain:
    def __init__(self, working_dir, trees, log_file=None, summary=None, name: str = "MC"):
        if isinstance(name, str):
            self.name = name
        else:
            raise ValueError("Unrecognized type for name!")
        # Setting the Working directory
        if isinstance(working_dir, str):
            if os.path.isdir(working_dir):
                self.working_dir = working_dir
            else:
                if isinstance(working_dir, str):
                    raise NotADirectoryError(working_dir)

                raise ValueError(working_dir)
        else:
            raise ValueError("Working directory type not recognized!")

        for new_dir in ["data", "plots", "clustering"]:
            try:
                os.mkdir(f"{self.working_dir}/{new_dir}")
            except FileExistsError:
                pass

        # Setting trees
        if isinstance(trees, TimeTreeSet):
            self.trees = trees
        elif isinstance(trees, str):
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
            if isinstance(log_file, str):
                raise FileNotFoundError(log_file)
            if log_file is not None:
                raise ValueError(log_file)

        if summary is not None:
            # Setting the summary tree
            if isinstance(summary, TimeTreeSet):
                self.summary = summary
            elif isinstance(summary, str):
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
        _matrix_file_name = (f"{self.working_dir}/data/"
                             f"{self.name if name == '' else name}"
                             f"{'_rf' if rf else ''}")
        if not os.path.exists(f"{_matrix_file_name}.{'csv.gz' if csv else 'npy'}"):
            dm = calc_pw_distances(self.trees, rf=rf)
            if csv:
                np.savetxt(fname=f"{_matrix_file_name}.csv.gz",
                           X=dm, delimiter=',', fmt='%i')
            else:
                np.save(file=f"{_matrix_file_name}.npy", arr=dm)
            return dm

        if csv:
            return np.genfromtxt(fname=f"{_matrix_file_name}.csv.gz",
                                 delimiter=',', dtype=int)

        return np.load(file=f"{_matrix_file_name}.npy")

    def get_key_names(self):
        """
        :return: list of all log file entry keys
        """
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
            return calc_ess(
                list(trace=self.log_data[ess_key][lower_i:(upper_i + 1)].dropna()))
        else:
            raise ValueError("Not (yet) implemented!")

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

        # todo here we can use DIST Literal...
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

    def split_trees_from_clustering(self, k, _overwrite=False):

        # todo change to the get_clustering funciton
        #  and also include the type of clustering that is wanted here
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

    @validate_literal_args(clustering_type=CLUSTERING_TYPE, dist_type=DIST)
    def get_clustree(self, k,
                     cluster_type: CLUSTERING_TYPE = "spectral",
                     dist_type: DIST = "rnni",
                     _overwrite: bool = False,
                     **kwargs):

        if not isinstance(k, int):
            raise TypeError("k must be an integer!")
        if k < 1:
            raise ValueError(f"Value {k} not supported!")

        if k == 1:
            return np.zeros(len(self.trees), dtype=int)
        if k > 1:
            match cluster_type:
                case "spectral":
                    beta = kwargs.get("beta", 1)
                    if not isinstance(beta, (float, int)):
                        raise TypeError("beta must be a float or int!")
                    _cluster_file_id = (f"{cluster_type}"
                                        f"{'' if float(beta) == 1.0 else f'-b{beta:.4g}'}")

                    _cluster_file_name = (f"{os.path.join(self.working_dir, 'data')}/"
                                          f"{self.name}_{_cluster_file_id}"
                                          f"-C{k}"
                                          f"-{dist_type}")

                    # WIP: this will be repetitive with more methods of clustering, put to fnct
                    if _overwrite:
                        try:
                            os.remove(f"{_cluster_file_name}.npy")
                        except FileNotFoundError:
                            pass

                    if os.path.exists(
                            f"{_cluster_file_name}.npy"):
                        clustering = np.load(f"{_cluster_file_name}.npy")
                    else:
                        # todo here we need to pass dist_type on to pwd_matrix
                        warnings.warn("Disttype is currently not implemented! WIP",
                                      stacklevel=3)
                        clustering = spectral_clustree_dm(self.pwd_matrix(),
                                                          n_clus=k, beta=beta)
                        np.save(file=_cluster_file_name, arr=clustering)

                    return clustering

                case _:
                    raise NotImplementedError(f"Unrecognized cluster type...{cluster_type}!")

    @validate_literal_args(mds_type=MDS_TYPES, dist_type=DIST)
    def get_mds_coords(self, mds_type: MDS_TYPES = 'tsne',
                       dim: int = 2, dist_type: DIST = 'rnni',
                       _overwrite: bool = False) -> np.ndarray:

        warnings.warn("Currently not fully implemented.. "
                      "Will only compute RNNI TSNE 2dim MDS.",
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

    def evaluate_clustering(self, kind="bic", _overwrite_clustering=False, add_random=True,
                            local=False):

        # todo parameters for overwrite currently not used

        # todo kind should be an enum
        if kind not in ["silhouette", "bic"]:
            raise ValueError(f"Kind {kind} not supported, choose either 'Silhouette' or 'BIC'.")

        # todo if _overwrite_plot delete plot and silhouette score
        # todo if _recalculate_clsutering delete clustering
        #  and also the plot and save file for the scores

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

            # _overwrite=True will recalculate the clustering too,
            # as this is being passed down to the bic() function
            plot_bic(treeset=self.trees,
                     clusterings=clusterings,
                     summaries_dict=summaries_dict,
                     max_cluster=max_cluster,
                     local_norm=local_norm,
                     add_random=add_random,
                     file_name=os.path.join(self.working_dir, "plots",
                                            f"BIC_{'random_' if add_random else ''}"
                                            f"{'local_' if local_norm else ''}{self.name}.pdf"))

    def plot_mds(self, mds_type: MDS_TYPES = 'tsne', dim: int = 2,
                 dist_type: DIST = 'rnni', plot_options: PlotOptions = None) -> None:
        """
        (WIP) Plotting simple MDS of a chain without any clustering.

        :param plot_options:
        :param mds_type: Type of MDS to plot
        :param dim: Dimension to plot to (only dim=2 supported atm)
        :param dist_type: Type of tree distance to use (RNNI or RF)
        :return: None
        """

        if dim != 2:
            raise ValueError("Unsupported dimension for MDS -- only supports dim=2 at the moment.")

        coords = self.get_mds_coords(mds_type=mds_type, dim=dim)
        if plot_options is None:
            plot_options = PlotOptions()
        self._fill_plot_defaults(plot_options, mds_type=mds_type, dist_type=dist_type)

        plot_coords(coords=coords, dim=dim, options=plot_options)

    def _fill_plot_defaults(self,
                            plot_options: PlotOptions,
                            mds_type: MDS_TYPES,
                            dist_type: DIST):
        """
        Filling default plotting options using plot_options.

        :param plot_options: The plot options to fill in
        :param mds_type: Type of MDS
        :param dist_type: Distance type
        :return: None
        """
        if not plot_options.was_explicit("filename"):
            plot_options.filename = os.path.join(self.working_dir, "plots",
                                                 f"{self.name}_{mds_type}-{dist_type}.pdf")
        if not plot_options.was_explicit("colors"):
            plot_options.colors = np.ones(len(self), dtype=int)
        if not plot_options.was_explicit("label"):
            plot_options.label = self.name
        if not plot_options.was_explicit("title"):
            plot_options.title = f"{self.name} - {mds_type}-{dist_type} MDS Plot"
