

from . import _variations
from ._constants import VAR_LIST
from treeoclock.trees.time_trees import TimeTreeSet


# TODO Tests and Documentation


class Centroid:
    """
    Class to compute the centroid algorithm with a set of parameters
    """

    def __init__(self, variation="inc-sub"):
        self.variation = variation


    def compute_centroid(self, trees: TimeTreeSet):

        # TODO this error may move to the actual centroid function
        if self.variation not in VAR_LIST:
            raise ValueError(f"The 'variation' parameter should be"
                             f"in {VAR_LIST} but {self.variation} was given.")

        return getattr(_variations, self.variation)()

        return 0
