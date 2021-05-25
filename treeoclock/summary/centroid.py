

from . import _variations
from ._constants import VAR_LIST
from treeoclock.trees.time_trees import TimeTreeSet


# TODO Tests and Documentation


class Centroid:
    """
    Class to compute the centroid algorithm with a set of parameters
    """

    def __init__(self, variation="inc-sub", n_cores=None):
        self.variation = variation
        self.n_cores = n_cores

    def compute_centroid(self, trees: TimeTreeSet):

        # TODO this error may move to the actual centroid function
        if self.variation not in VAR_LIST:
            raise ValueError(f"The 'variation' parameter should be"
                             f" in {VAR_LIST} but {self.variation} was given.")
        if not isinstance(self.n_cores, (int)):
            raise ValueError(f"The 'n_cores' parameter should be an integer"
                             f" but {self.n_cores} was given.")

        # TODO the starting  tree needs a parameter list [fm, efm, last, random,]

        # TODO Argument trees need to be passed
        return getattr(_variations, self.variation)(trees)

        return 0
