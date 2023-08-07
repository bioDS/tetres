import os
import numpy as np

from tetres.judgment.chain import Chain
from tetres.judgment._pairwise_distance_matrix import calc_pw_distances_two_sets
from tetres.judgment import _gelman_rubin_diag as grd
from tetres.judgment import _cladesetcomparator as csc
from tetres.summary.centroid import Centroid
from tetres.summary.annotate_centroid import annotate_centroid
from tetres.judgment._discrete_cladesetcomparator import discrete_cladeset_comparator
from tetres.judgment._extract_cutoff import _extract_cutoff
from tetres.judgment.burnin_detection import burn_detector


class MultiChain():
    def __init__(self, m_chains, trees, log_files, working_dir, name="cMC"):
        # todo currently still missing the summary parameter of MChain, but its not really used at this point
        self.m_chains = m_chains
        if self.m_chains < 2 or self.m_chains < 0:
            raise ValueError("Wrong usage of MultiChain")
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
            if len(trees) != self.m_chains or len(log_files) != self.m_chains:
                raise ValueError("m_chains and number of trees do not match!")
            self.tree_files = trees.copy()
            self.log_files = log_files.copy()
            for i in range(self.m_chains):
                self.MChain_list.append(Chain(trees=trees[i], log_file=log_files[i],
                                              working_dir=working_dir, name=f"{self.name}_{i}"))

        elif type(trees) is str:
            if not os.path.exists(f"{self.working_dir}/{trees}"):
                raise FileNotFoundError(f"Given trees file {self.working_dir}/{trees} does not exist!")
            if not os.path.exists(f"{self.working_dir}/{log_files}"):
                raise FileNotFoundError(f"Given trees file {self.working_dir}/{log_files} does not exist!")
            self.MChain_list.append(
                Chain(trees=trees, log_file=log_files, working_dir=working_dir, name=f"{self.name}_{0}"))
            for i in range(1, self.m_chains):
                self.MChain_list.append(
                    Chain(trees=f"chain{i}{trees}", log_file=f"chain{i}{log_files}",
                          working_dir=working_dir, name=f"{self.name}_{i}"))
        else:
            raise ValueError("Unrecognized argument types of trees and log_files!")
        if not all([self.MChain_list[0].trees.map == self.MChain_list[c].trees.map for c in range(1, self.m_chains)]):
            raise ValueError(
                "The taxon maps in the tree files are not the same, this can lead to problems! Fix before continue!")

    def pwd_matrix(self, index1, index2=None, csv: bool = False, rf: bool = False):
        if type(index1) is not int:
            raise ValueError("Unrecognized index type!")
        if index1 >= self.m_chains:
            raise IndexError("Given Index out of range!")

        if index2 is None:
            return self.MChain_list[index1].pwd_matrix(name=self.MChain_list[index1].name, csv=csv, rf=rf)
        else:
            if type(index2) is not int:
                raise ValueError("Unrecognized index type!")
            if index1 == index2:
                raise ValueError("Indeces are equal, only call with one index for pairwise distances of one set!")
            if index2 >= self.m_chains:
                raise IndexError("Given Index out of range!")
            if index2 < index1:
                index1, index2 = index2, index1
            if not os.path.exists(
                    f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.{'csv.gz' if csv else 'npy'}"):
                dm = calc_pw_distances_two_sets(self.MChain_list[index1].trees, self.MChain_list[index2].trees, rf=rf)
                if csv:
                    np.savetxt(
                        fname=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.csv.gz",
                        X=dm, delimiter=',',
                        fmt='%i')
                else:
                    np.save(file=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.npy",
                            arr=dm)
                return dm
            else:
                if csv:
                    return np.genfromtxt(
                        fname=f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.csv.gz",
                        delimiter=',', dtype=int)
                else:
                    return np.load(f"{self.working_dir}/data/{self.name}_{index1}_{index2}{'_rf' if rf else ''}.npy")

    def pwd_matrix_all(self, rf: bool = False):
        if not os.path.exists(f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy"):
            combined_matrix = np.array([])
            for i in range(self.m_chains):
                cur_row = np.array([])
                for j in range(self.m_chains):
                    # print(i, j)
                    if i < j:
                        cur_row = np.concatenate((cur_row, self.pwd_matrix(i, j, rf=rf)),
                                                 axis=1) if cur_row.size else self.pwd_matrix(i, j, rf=rf)
                    elif i > j:
                        cur_row = np.concatenate((cur_row, np.zeros((len(self[i].trees), len(self[j].trees)))),
                                                 axis=1) if cur_row.size else np.zeros(
                            (len(self[i].trees), len(self[j].trees)))
                    elif i == j:
                        cur_row = np.concatenate((cur_row, self.pwd_matrix(i, rf=rf)),
                                                 axis=1) if cur_row.size else self.pwd_matrix(i, rf=rf)
                    # print(cur_row.shape)
                combined_matrix = np.concatenate((combined_matrix, cur_row),
                                                 axis=0) if combined_matrix.size else cur_row
            np.save(file=f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy", arr=combined_matrix)
            return combined_matrix
        else:
            return np.load(f"{self.working_dir}/data/{self.name}_all{'_rf' if rf else ''}.npy")

    def extract_cutoff(self, i, j, ess_threshold=200, smoothing=0.6, gr_boundary=0.02, smoothing_average="mean",
                       subsampling=False, _overwrite=False, burnin=0):
        # The tree and log file strings for storing the cutoff
        tree_file = f"{self.working_dir}/cutoff_files/{self[i].name}_{self[j].name}{f'_burnin-{burnin}' if burnin != 0 else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}{f'_subsample-{subsampling}' if subsampling else ''}_smoothing-{smoothing}_{smoothing_average}_boundary-{gr_boundary}.trees"
        log_file = f"{self.working_dir}/cutoff_files/{self[i].name}_{self[j].name}{f'_burnin-{burnin}' if burnin != 0 else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}{f'_subsample-{subsampling}' if subsampling else ''}_smoothing-{smoothing}_{smoothing_average}_boundary-{gr_boundary}.log"
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
        start, end = self.gelman_rubin_cut(i=i, j=j, smoothing=smoothing, ess_threshold=ess_threshold,
                                           smoothing_average=smoothing_average, _gr_boundary=gr_boundary,
                                           _subsampling=subsampling, _overwrite=_overwrite, burnin=burnin)
        # Check if no cutoff raise ValueError
        if start < 0 or end < 1:
            raise ValueError("No cutoff exist, therefore no cutoff can be extracted!")

        _extract_cutoff(multichain=self, i=i, start=start, end=end, tree_file=tree_file, log_file=log_file)

    def gelman_rubin_all_chains_density_plot(self, samples: int = 100):
        return grd.gelman_rubin_all_chains_density_plot(self, samples=samples)

    def psrf_density_trace_plot(self, interval, i=0, j=1):
        return grd.density_trace_plot(self, interval, i=i, j=j)

    def gelman_rubin_parameter_choice_plot(self, i, j, _subsampling=False, _gr_boundary=0.02, smoothing_average="mean"):
        return grd.gelman_rubin_parameter_choice_plot(self, i, j, _subsampling=_subsampling, _gr_boundary=_gr_boundary,
                                                      smoothing_average=smoothing_average)

    def gelman_rubin_cut(self, i, j, smoothing=0.5, ess_threshold=200, pseudo_ess_range=100, _overwrite=False,
                         smoothing_average="mean", _subsampling=False, _gr_boundary=0.02, burnin=0):

        # sorting to make sense for dm_ij interpretation
        if j < i:
            i, j = j, i

        # First if overwrite delete previous computation file
        if _overwrite:
            try:
                os.remove(
                    f"{self.working_dir}/data/{self.name}_{i}_{j}_gelman_rubin_cutoff{f'_burnin_{burnin}' if burnin != 0 else ''}{f'_subsampling-{_subsampling}' if _subsampling else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}_smoothing-{smoothing}_{smoothing_average}_boundary-{_gr_boundary}")
            except FileNotFoundError:
                pass
        # IF no overwrite and the file exists simply read the values from the file
        if os.path.exists(
                f"{self.working_dir}/data/{self.name}_{i}_{j}_gelman_rubin_cutoff{f'_burnin_{burnin}' if burnin != 0 else ''}{f'_subsampling-{_subsampling}' if _subsampling else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}_smoothing-{smoothing}_{smoothing_average}_boundary-{_gr_boundary}"):
            with open(
                    f"{self.working_dir}/data/{self.name}_{i}_{j}_gelman_rubin_cutoff{f'_burnin_{burnin}' if burnin != 0 else ''}{f'_subsampling-{_subsampling}' if _subsampling else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}_smoothing-{smoothing}_{smoothing_average}_boundary-{_gr_boundary}",
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
                                                  _gr_boundary=_gr_boundary,
                                                  burnin=burnin)
        # Write the cutoff boundaries to a file, if it already exists skip this part
        if _subsampling:
            # If subsampling the start and end need to be adapted, each subsampling is dividing the chain length by 2
            for _ in range(_subsampling):
                cut_start, cut_end = cut_start * 2, cut_end * 2
        try:
            # Write the cut start and end to a file, if it exists raises an Error
            with open(
                    f"{self.working_dir}/data/{self.name}_{i}_{j}_gelman_rubin_cutoff{f'_burnin_{burnin}' if burnin != 0 else ''}{f'_subsampling-{_subsampling}' if _subsampling else ''}{'' if ess_threshold == 0 else f'_ess-{ess_threshold}'}_smoothing-{smoothing}_{smoothing_average}_boundary-{_gr_boundary}",
                    "x") as f:
                f.write(f"{cut_start}\n{cut_end}")
        except FileExistsError:
            raise Exception(
                "The file for writing the GRD cut start and end exists already, this shouldn't happen here!")
        # return the start and end values
        return cut_start, cut_end

    def detect_burnin(self, i, j, traces=["posterior", "likelihood", "prior"]):
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
        return discrete_cladeset_comparator(tree_set_i=self[i].trees, tree_set_j=self[j].trees, plot=plot,
                                            burnin=burnin,
                                            file=f"{self.working_dir}/plots/{self.name}_discrete_cc_{i}_{j}_burn-{burnin}.png")

    def cen_for_each_chain(self, repeat=5):
        raise ValueError("Current Work in Progress and not supported for now.")
        # todo should be its own script, adhering to the naming conventions to fit in with clustering
        #  clustering stuff should use the output from this funciton for 1 clustering
        for _ in range(repeat):
            for chain in self.MChain_list:
                cen = Centroid()
                cen_tree, sos = cen.compute_centroid(chain.trees)
                cen_tree = annotate_centroid(cen_tree, chain.trees)
                try:
                    file = open(f"{chain.working_dir}/data/{chain.name}_cen_sos.log", "r")
                    cur_sos = int(file.readline())
                    file.close()
                except FileNotFoundError:
                    # Setting cur_sos so that the folowing if will be Treu
                    cur_sos = sos + 1
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
                    cen_tree.write_nexus(chain.trees.map, file_name=f"{chain.working_dir}/{chain.name}_cen.tree",
                                         name=f"Cen_{chain.name}")
