import os
import warnings

import numpy as np

from tetres.clustree.spectral_clustree import spectral_clustree_dm
from tetres.judgement import _cladesetcomparator as csc
from tetres.judgement import _gelman_rubin_diag as grd
from tetres.judgement._discrete_cladesetcomparator import discrete_cladeset_comparator
from tetres.judgement._extract_cutoff import _extract_cutoff
from tetres.judgement._pairwise_distance_matrix import calc_pw_distances_two_sets
from tetres.judgement.burnin_detection import burn_detector
from tetres.judgement.chain import Chain
from tetres.utils.decorators import validate_literal_args
from tetres.utils.literals import DIST, MDS_TYPES, CLUSTERING_TYPE
from tetres.visualize.mds_coord_compuation import _tsne_coords_from_pwd
from tetres.visualize.plot_config import PlotOptions
from tetres.visualize.plot_coords import plot_coords


class MultiChain():
    def __init__(self, m_chains, trees, log_files, working_dir, name="cMC"):
        self.m_chains = m_chains
        if self.m_chains < 2 or self.m_chains < 0:
            raise ValueError("Wrong usage of MultiChain")
        self.name = name

        # todo was unreachable at the end of the if statements...
        # if not all(
        #         [self.MChain_list[0].trees.map == self.MChain_list[c].trees.map
        #          for c in range(1, self.m_chains)]):
        #     raise ValueError(
        #         "The taxon maps in the tree files are not the same,"
        #         " this can lead to problems! Fix before continue!")

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
            if len(trees) != self.m_chains or len(log_files) != self.m_chains:
                raise ValueError("m_chains and number of trees do not match!")
            self.tree_files = trees.copy()
            self.log_files = log_files.copy()
            for i in range(self.m_chains):
                self.MChain_list.append(Chain(trees=trees[i], log_file=log_files[i],
                                              working_dir=working_dir,
                                              name=f"{self.name}_chain{i}"))

        elif type(trees) is str:
            if not os.path.exists(f"{self.working_dir}/{trees}"):
                raise FileNotFoundError(
                    f"Given trees file {self.working_dir}/{trees} does not exist!")
            if not os.path.exists(f"{self.working_dir}/{log_files}"):
                raise FileNotFoundError(
                    f"Given trees file {self.working_dir}/{log_files} does not exist!")
            self.MChain_list.append(
                Chain(trees=trees, log_file=log_files, working_dir=working_dir,
                      name=f"{self.name}_chain{0}"))
            for i in range(1, self.m_chains):
                self.MChain_list.append(
                    Chain(trees=f"chain{i}{trees}", log_file=f"chain{i}{log_files}",
                          working_dir=working_dir,
                          name=f"{self.name}_chain{i}"))
        else:
            raise ValueError("Unrecognized argument types of trees and log_files!")

    def pwd_matrix(self, index1, index2=None, csv: bool = False, rf: bool = False):
        # todo make this use a target and return the "all" matrix if requested to avoid duplicate
        #  name
        if type(index1) is not int:
            raise ValueError("Unrecognized index type!")
        if index1 >= self.m_chains:
            raise IndexError("Given Index out of range!")

        if index2 is None:
            return self.MChain_list[index1].pwd_matrix(name=self.MChain_list[index1].name, csv=csv,
                                                       rf=rf)
        else:
            if type(index2) is not int:
                raise ValueError("Unrecognized index type!")
            if index1 == index2:
                raise ValueError(
                    "Indeces are equal, only call with one index for pairwise distances of one set!")
            if index2 >= self.m_chains:
                raise IndexError("Given Index out of range!")
            if index2 < index1:
                index1, index2 = index2, index1
            if not os.path.exists(
                    f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.{'csv.gz' if csv else 'npy'}"):
                dm = calc_pw_distances_two_sets(self.MChain_list[index1].trees,
                                                self.MChain_list[index2].trees, rf=rf)
                if csv:
                    np.savetxt(
                        fname=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.csv.gz",
                        X=dm, delimiter=',',
                        fmt='%i')
                else:
                    np.save(
                        file=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.npy",
                        arr=dm)
                return dm
            else:
                if csv:
                    return np.genfromtxt(
                        fname=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.csv.gz",
                        delimiter=',', dtype=int)
                else:
                    return np.load(
                        f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.npy")

    def pwd_matrix_all(self, rf: bool = False):
        if not os.path.exists(f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy"):
            combined_matrix = np.array([])
            for i in range(self.m_chains):
                cur_row = np.array([])
                for j in range(self.m_chains):
                    # print(i, j)
                    if i < j:
                        cur_row = np.concatenate((cur_row, self.pwd_matrix(i, j, rf=rf)),
                                                 axis=1) if cur_row.size else self.pwd_matrix(i, j,
                                                                                              rf=rf)
                    elif i > j:
                        cur_row = np.concatenate(
                            (cur_row, np.zeros((len(self[i].trees), len(self[j].trees)))),
                            axis=1) if cur_row.size else np.zeros(
                            (len(self[i].trees), len(self[j].trees)))
                    elif i == j:
                        cur_row = np.concatenate((cur_row, self.pwd_matrix(i, rf=rf)),
                                                 axis=1) if cur_row.size else self.pwd_matrix(i,
                                                                                              rf=rf)
                    # print(cur_row.shape)
                combined_matrix = np.concatenate((combined_matrix, cur_row),
                                                 axis=0) if combined_matrix.size else cur_row
            np.save(file=f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy",
                    arr=combined_matrix)
            return combined_matrix
        else:
            return np.load(f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy")

    def extract_cutoff(self, i, j, ess_threshold=200, smoothing=1,
                       tolerance=0.02, smoothing_average="mean",
                       subsampling=False, _overwrite=False, burnin=0):
        # The tree and log file strings for storing the cutoff
        _base_file = (f"{self.working_dir}/cutoff_files/"
                      f"{self[i].name}_{self[j].name}"
                      f"{f'_burnin-{burnin}' if burnin != 0 else ''}"
                      f"{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}"
                      f"{f'_subsample-{subsampling}' if subsampling else ''}"
                      f"_smoothing-{smoothing}_{smoothing_average}_boundary-{tolerance}")
        tree_file = f"{_base_file}.trees"
        log_file = f"{_base_file}.log"
        if os.path.exists(tree_file) and os.path.exists(log_file) and not _overwrite:
            # Files have already been extracted, nothing will be done
            return 0
        else:
            try:
                os.remove(tree_file)
            except FileNotFoundError:
                pass
            try:
                os.remove(log_file)
            except FileNotFoundError:
                pass
        # Try creating the cutoff directory
        try:
            os.mkdir(f"{self.working_dir}/cutoff_files")
        except FileExistsError:
            pass
        # computing start and end from the current Multichain
        start, end = self.gelman_rubin_cut(i=i, j=j, smoothing=smoothing,
                                           ess_threshold=ess_threshold,
                                           smoothing_average=smoothing_average, tolerance=tolerance,
                                           _subsampling=subsampling, _overwrite=_overwrite,
                                           burnin=burnin)
        # Check if no cutoff raise ValueError
        if start < 0 or end < 1:
            raise ValueError("No cutoff exist, therefore no cutoff can be extracted!")

        _extract_cutoff(multichain=self, i=i, start=start, end=end, tree_file=tree_file,
                        log_file=log_file)

    def gelman_rubin_all_chains_density_plot(self, samples: int = 100):
        return grd.gelman_rubin_all_chains_density_plot(self, samples=samples)

    def psrf_density_trace_plot(self, interval, i=0, j=1, no_smooth=False):
        return grd.density_trace_plot(self, interval, i=i, j=j, no_smooth=no_smooth)

    def gelman_rubin_parameter_choice_plot(self, i, j, _subsampling=False, tolerance=0.02,
                                           smoothing_average="mean"):
        return grd.gelman_rubin_parameter_choice_plot(self, i, j, _subsampling=_subsampling,
                                                      tolerance=tolerance,
                                                      smoothing_average=smoothing_average)

    def gelman_rubin_cut(self, i, j, smoothing=1, ess_threshold=200, pseudo_ess_range=100,
                         _overwrite=False, smoothing_average="mean", _subsampling=False,
                         tolerance=0.02, burnin=0):

        # sorting to make sense for dm_ij interpretation
        if j < i:
            i, j = j, i

        cutoff_file_name = (f"{self.working_dir}/"
                            f"data/"
                            f"{self.name}_{i}_{j}_"
                            f"gelman_rubin_cutoff"
                            f"{f'_burnin_{burnin}' if burnin != 0 else ''}"
                            f"{f'_subsampling-{_subsampling}' if _subsampling else ''}"
                            f"{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}"
                            f"_smoothing-{smoothing}_{smoothing_average}_boundary-{tolerance}")

        # First if overwrite delete previous computation file
        if _overwrite:
            try:
                os.remove(cutoff_file_name)
            except FileNotFoundError:
                pass
        # IF no overwrite and the file exists simply read the values from the file
        if os.path.exists(cutoff_file_name):
            with open(
                    cutoff_file_name,
                    "r") as file:
                cut_start = int(file.readline())
                cut_end = int(file.readline())
            return cut_start, cut_end
        # File did not exist or was deleted
        # Computing start and end of the cut via the GRD function
        cut_start, cut_end = grd.gelman_rubin_cut(self, i=i, j=j,
                                                  smoothing=smoothing,
                                                  ess_threshold=ess_threshold,
                                                  pseudo_ess_range=pseudo_ess_range,
                                                  smoothing_average=smoothing_average,
                                                  _subsampling=_subsampling,
                                                  tolerance=tolerance,
                                                  burnin=burnin)
        # Write the cutoff boundaries to a file, if it already exists skip this part
        if _subsampling:
            # If subsampling the start and end need to be adapted,
            # each subsampling is dividing the chain length by 2
            for _ in range(_subsampling):
                cut_start, cut_end = cut_start * 2, cut_end * 2
        try:
            # Write the cut start and end to a file, if it exists raises an Error
            with open(cutoff_file_name, "x") as f:
                f.write(f"{cut_start}\n{cut_end}")
        except FileExistsError:
            raise Exception(
                "The file for writing the GRD cut start and end exists already. "
                "This shouldn't happen here!")
        # return the start and end values
        return cut_start, cut_end

    def detect_burnin(self, i, j, traces=None):
        if traces is None:
            traces = ["posterior", "likelihood", "prior"]
        return burn_detector(self, chain_i=i, chain_j=j, traces=traces)

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
        # Calling the BEAST applauncher cladesetcomparator function with a wrap for more than 2 chains
        csc.cladesetcomp(self, beast_applauncher, burnin=burnin)

    def discrete_cladesetcomparator(self, i, j, plot=True, burnin=0):
        # Calling the discrete version of the cladesetcomparator, return the single value
        return discrete_cladeset_comparator(tree_set_i=self[i].trees, tree_set_j=self[j].trees,
                                            plot=plot,
                                            burnin=burnin,
                                            file=f"{self.working_dir}/plots/{self.name}_discrete_cc_{i}_{j}_burn-{burnin}.png")

    @validate_literal_args(mds_type=MDS_TYPES, dist_type=DIST)
    def get_mds_coords(self, target="all", mds_type: MDS_TYPES = 'tsne',
                       dim: int = 2, dist_type: DIST = 'rnni',
                       _overwrite: bool = False, seed: (None, int) = None) -> np.ndarray:
        warnings.warn("Currently not fully implemented.. Will only compute RNNI TSNE 2dim MDS.",
                      stacklevel=3)

        # todo will need to add more possible parameters, kwargs... beta for spectral

        match target:
            case int() as i:
                return self[i].get_mds_coords(mds_type=mds_type, dim=dim,
                                              dist_type=dist_type, _overwrite=_overwrite,
                                              seed=seed)
            case "all":
                mds_coords_filename = os.path.join(self.working_dir, "data",
                                                   f"{self.name}_{mds_type}-{dist_type}"
                                                   f"{"" if seed is None else "_seed-" + str(seed)}"
                                                   f".npy")

                if not os.path.exists(mds_coords_filename) or _overwrite:
                    # need to first calculate the coordinates with the correct function

                    # todo no really using the mds_type atm, need to implement that
                    coords = _tsne_coords_from_pwd(self.pwd_matrix_all(rf=False), dim=dim,
                                                   seed=seed)
                    np.save(file=mds_coords_filename, arr=coords)
                else:
                    # have to read the already computed coords
                    coords = np.load(file=mds_coords_filename)
                return coords
            case _:
                raise ValueError(f"Invalid target type '{target}'. "
                                 f"Choose from: 'all, 0 < index < {len(self.MChain_list)}'")

    def plot_mds(self, target="all", mds_type: MDS_TYPES = 'tsne',
                 dim: int = 2, dist_type: DIST = 'rnni',
                 plot_options: PlotOptions = None, seed: (None, int) = None,
                 _overwrite: bool = False) -> None:
        """
        (WIP) Plotting simple MDS of MutliChain, either single chain or all together.

        :param _overwrite: If true recompute the MDS coordinates
        :param plot_options: Plot options, see plot_config.py for more infos
        :param seed: Random number seed for tSNE embedding
        :param target: Either 'all' or index for which chain to plot
        :param mds_type: Type of MDS coords to plot
        :param dim: Dimension to plot to (only dim=2 supported atm)
        :param dist_type: Type of tree distance to use (RNNI or RF)
        :return: None
        """

        match target:
            case int() as i:
                return self[i].plot_mds(mds_type=mds_type, dim=dim, dist_type=dist_type,
                                        _overwrite=_overwrite, seed=seed)
            case "all":
                if dim != 2:
                    raise ValueError("Currently not supported dimension, only dim=2 accepted.")

                if plot_options is None:
                    plot_options = PlotOptions()
                self._fill_plot_defaults(plot_options, mds_type=mds_type, dist_type=dist_type)

                coords = self.get_mds_coords(target=target, mds_type=mds_type, dim=dim, seed=seed,
                                             _overwrite=_overwrite)

                plot_coords(coords=coords, dim=dim, options=plot_options)
                return None

            case _:
                raise ValueError(f"Invalid target type '{target}'. "
                                 f"Choose from: 'all, 0 < index < {len(self.MChain_list)}'")

    def _fill_plot_defaults(self,
                            plot_options: PlotOptions,
                            mds_type: MDS_TYPES,
                            dist_type: DIST):
        """
        Filling in default plotting options using plot_options.

        :param plot_options: The plot options to fill in
        :param mds_type: Type of MDS
        :param dist_type: Distance type
        :return: None
        """
        if not plot_options.was_explicit("filename"):
            plot_options.filename = os.path.join(self.working_dir, "plots",
                                                 f"{self.name}_{mds_type}-{dist_type}.pdf")
        if not plot_options.was_explicit("colors"):
            plot_options.colors = np.concatenate(
                [[i] * len(chain) for i, chain in enumerate(self.MChain_list)])
        if not plot_options.was_explicit("label"):
            plot_options.label = {i: c.name for i, c in enumerate(self.MChain_list)}
        if not plot_options.was_explicit("title"):
            plot_options.title = f"{self.name} - {mds_type}-{dist_type} MDS Plot"

    @validate_literal_args(mds_type=MDS_TYPES, dist_type=DIST)
    def get_clustree(self, target="all", k: int = 1,
                     cluster_type: CLUSTERING_TYPE = "spectral",
                     dist_type: DIST = "rnni",
                     _overwrite: bool = False,
                     **kwargs):

        if not isinstance(k, int):
            raise TypeError("k must be an integer!")
        if k < 1:
            raise ValueError(f"Value {k} not supported!")

        match target:
            case int() as i:
                return self[i].get_clustree(k=k,
                                            cluster_type=cluster_type,
                                            dist_type=dist_type,
                                            _overwrite=_overwrite,
                                            **kwargs)
            case "all":
                if k == 1:
                    total_len = sum(len(i) for i in self.MChain_list)
                    return np.zeros(total_len, dtype=int)
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

                            # WIP: this will be repetitive with more methods of clustering
                            #  --> put to fnct
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
                                warnings.warn("Dist_type is currently not implemented! WIP",
                                              stacklevel=3)
                                clustering = spectral_clustree_dm(self.pwd_matrix_all(),
                                                                  n_clus=k, beta=beta)
                                np.save(file=_cluster_file_name, arr=clustering)

                            return clustering

                        case _:
                            raise NotImplementedError(
                                f"Cluster Type {cluster_type} not implemented.")

            case _:
                raise ValueError(f"Invalid target type '{target}'.")
