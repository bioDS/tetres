# todo this is WIP with old funcitons that should be in this subpackage


# todo wip
    def plot_chains_with_summaries(self, rf: bool = False, dim: int = 2):
        new_cmchain = coupled_MChains(self.m_MChains, self.tree_files, self.log_files, self.working_dir, name=f"{self.name}_cen")
        new_cmchain.m_MChains += new_cmchain.m_MChains
        for chain in self.MChain_list:
            cur_chain = MChain(trees=TimeTreeSet(f"{chain.working_dir}/{chain.name}_cen.tree"), log_file=None,
                               working_dir=self.working_dir, name=f"{chain.name}_cen")
            new_cmchain.MChain_list.append(cur_chain)
        all_chains_added_summaries(new_cmchain, rf=rf, dim=dim)


def plot_clustree_all(self, n_clus=2, beta=1, rf: bool = False, dim: int = 2):
    all_chains_spectral_clustree(self, beta=beta, n_clus=n_clus, rf=rf, dim=dim)


def split_all_trees(self, n_clus, beta=1, burn_in=5):
    clustering = self.clustree_all(n_clus=n_clus, beta=beta)
    try:
        os.mkdir(f"{self.working_dir}/{n_clus}_cluster")  # todo maybe name the cluster based on rf also!
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
            for taxa in range(1, self[0].trees[0].ctree.num_leaves + 1):
                f.write(f"\t\t\t{self[0].trees.map[taxa]}\n")
            f.write("\t\t\t;\nEnd;\nBegin trees;\n\tTranslate\n")
            for taxa in range(1, self[0].trees[0].ctree.num_leaves):
                f.write(f"\t\t\t{taxa} {self[0].trees.map[taxa]},\n")
            f.write(
                f"\t\t\t{self[0].trees[0].ctree.num_leaves} {self[0].trees.map[self[0].trees[0].ctree.num_leaves]}\n")
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

                # offset = 10 + (2 * self[0].trees[0].ctree.num_leaves)  # this approach does not work when there is a random empty line inserted, should be all the same file format but all are different ...
                file = f"{self.working_dir}/{self.tree_files[chain]}"
                offset = 0
                for line in open(file):
                    if re_tree.match(line):
                        break
                    offset += 1

                line = index + 1 + offset
                cur_tree = linecache.getline(file, line)
                cur_tree = f'{re.split(re_tree, cur_tree)[1][:re.split(re_tree, cur_tree)[1].rfind(")") + 1]};'

                cur_file.write(f"tree STATE_{sample_dict[clustering[index]]} = {cur_tree}\n")

                cur_file.close()
                # increasing the sample for the respective clustering
                sample_dict[clustering[index]] += 1
    for c in range(n_clus):
        cur_file = open(f"{self.working_dir}/{n_clus}_cluster/clus_{c}.trees", "a")
        cur_file.write("End;")
        cur_file.close()
    return 0


from tetres.clustree.spectral_clustree import _spectral_clustree
import numpy as np
import os
from tetres.visualize.plot_coords import plot_coords
from tetres.trees.time_trees import TimeTreeSet


# todo include the plot all chains tsne funciton here at some point!


def all_chains_spectral_clustree(cmchain, beta=1, n_clus=2, rf: bool = False, dim: int = 2):
    clustering = cmchain.clustree_all(n_clus=n_clus, beta=beta, rf=rf)
    coords = cmchain.tsne_all(rf=rf, dim=dim)
    plot_coords(coords, colors=clustering, filename=f"{cmchain.working_dir}/plots/{cmchain.name}_"
                                                    f"{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all"
                                                    f"{'_rf' if rf else ''}.png")


def all_chains_added_summaries(cmchain, rf: bool = False, dim: int = 2):
    coords = cmchain.tsne_all(rf=rf, dim=dim)

    clustering = np.ones(coords.shape[0])
    # Todo hard coded for one center per chain i.e. half the chains in the given cmchain are centers!
    for i in range(1, int(cmchain.m_MChains / 2) + 1):
        clustering[-i] = i + 1

    plot_coords(coords, colors=clustering, filename=f"{cmchain.working_dir}/plots/{cmchain.name}_"
                                                    f"summaries{'_rf' if rf else ''}.png",
                centers=int(cmchain.m_MChains / 2))


    def similarity_matrix_all(self, beta=1, rf: bool = False):
        if not os.path.exists(f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity_all{'_rf' if rf else ''}.npy"):
            matrix = self.pwd_matrix_all(rf=rf)
            similarity = np.exp(-beta * matrix / matrix.std())
            np.save(file=f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity_all{'_rf' if rf else ''}.npy",
                    arr=similarity)
            return similarity
        else:
            return np.load(f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}similarity_all{'_rf' if rf else ''}.npy")

    def clustree_all(self, n_clus=2, beta=1, rf: bool = False):
        if n_clus != 1:
            if not os.path.exists(
                    f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all{'_rf' if rf else ''}.npy"):
                similarity = self.similarity_matrix_all(beta=beta, rf=rf)
                similarity = similarity + similarity.transpose()
                clustering = _spectral_clustree(similarity, n_clus=n_clus)
                np.save(
                    file=f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all{'_rf' if rf else ''}.npy",
                    arr=clustering)
                return clustering
            else:
                return np.load(
                    f"{self.working_dir}/data/{self.name}_{'' if beta == 1 else f'{beta}_'}{n_clus}clustering_all{'_rf' if rf else ''}.npy")
        else:
            return np.ones(len(self[0].trees)*len(self))

    def tsne_all(self, rf: bool = False, dim: int = 2):
        # try:
        if os.path.exists(f"{self.working_dir}/data/{self.name}_{dim}D_tSNE_all{'_rf' if rf else ''}.npy"):
            coords = np.load(
                f"{self.working_dir}/data/{self.name}_{dim}D_tSNE_all{'_rf' if rf else ''}.npy")
        # except FileNotFoundError:
        else:
            coords, kl_divergence = _tsne_coords_from_pwd(pwd_matrix=self.pwd_matrix_all(rf=rf), dim=dim)
            np.save(
                file=f"{self.working_dir}/data/{self.name}_{dim}D_tSNE_all{'_rf' if rf else ''}.npy",
                arr=coords)
            # todo writing the kl divergence somewhere
        return coords