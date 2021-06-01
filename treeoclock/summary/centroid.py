from . import _variations
from ._constants import VAR_LIST, START_LIST
from treeoclock.trees.time_trees import TimeTreeSet


# TODO Tests and Documentation


class Centroid:
    """
    Class to compute the centroid algorithm with a set of parameters
    """

    def __init__(self, variation="inc-sub", n_cores=None, start='random'):
        self.variation = variation  # Which centroid variation to compute
        self.n_cores = n_cores  # How many cores to use whenever Multiprocessing is/will be used
        self.start = start  # Specifying which tree to choose in case of multiple options with the same quality
        # depending on the variation if it will be used

    def compute_centroid(self, trees: TimeTreeSet):

        # Checking if all parameters are correctly set up
        if self.variation not in VAR_LIST:
            raise ValueError(f"The 'variation' parameter should be"
                             f" in {VAR_LIST} but {self.variation} was given.")
        if not (isinstance(self.n_cores, int) or self.n_cores is None):
            raise ValueError(f"The 'n_cores' parameter should be an integer"
                             f" but {self.n_cores} was given.")
        if self.start not in START_LIST:
            raise ValueError(f"The 'variation' parameter should be"
                             f" in {START_LIST} but {self.start} was given.")

        return getattr(_variations, self.variation)(trees)

        return 0
