import _variations
from treeoclock.summary._constants import VAR_LIST, SELECT_LIST
from treeoclock.trees.time_trees import TimeTreeSet


# TODO Tests and Documentation


class Centroid:
    """
    Class to compute the centroid algorithm with a set of parameters
    """

    def __init__(self, variation="greedy", n_cores=None, select='random'):
        self.variation = variation  # Which centroid variation to compute
        self.n_cores = n_cores  # How many cores to use whenever Multiprocessing is/will be used
        self.select = select  # Specifying which tree to choose in case of multiple options with the same quality
        # depending on the variation if it will be used

    def compute_centroid(self, trees: TimeTreeSet):

        # Checking if all parameters are correctly set up
        if self.variation not in VAR_LIST:
            raise ValueError(f"The 'variation' parameter should be"
                             f" in {VAR_LIST} but {self.variation} was given.")
        if not (isinstance(self.n_cores, int) or self.n_cores is None):
            raise ValueError(f"The 'n_cores' parameter should be an integer"
                             f" but {self.n_cores} was given.")
        if self.select not in SELECT_LIST:
            raise ValueError(f"The 'select' parameter should be"
                             f" in {SELECT_LIST} but {self.select} was given.")

        return getattr(_variations, self.variation)(trees=trees, n_cores=self.n_cores, select=self.select)


if __name__ == '__main__':
    d_name = "Dengue"

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    mycen = Centroid()

    cen = mycen.compute_centroid(myts)

    print(cen)