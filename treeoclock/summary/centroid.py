from treeoclock.summary import _variations
from treeoclock.summary._constants import SELECT_LIST, START_LIST
from treeoclock.trees.time_trees import TimeTreeSet
from treeoclock.summary.frechet_mean import frechet_mean

import random

# TODO Tests and Documentation


class Centroid:
    """
    Class to compute the centroid algorithm with a set of parameters
    """

    def __init__(self, variation="greedy", n_cores=None, select='random', start='FM', subsample_size=200):
        self.variation = variation  # Which centroid variation to compute
        self.n_cores = n_cores  # How many cores to use whenever Multiprocessing is/will be used
        self.select = select  # Specifying which tree to choose in case of multiple options with the same quality
        self.start = start  # Specifying the starting position for the centroid algorithm
        self.subsample_size = subsample_size  # Size of the subsamples for some of the centroid variations

    def compute_centroid(self, trees: TimeTreeSet):

        # Checking if all parameters are correctly set up
        if not hasattr(_variations, self.variation):
            raise ValueError(f"The 'variation' parameter should be in _variations.py"
                             f" but {self.variation} was given.")
        if not (isinstance(self.n_cores, int) or self.n_cores is None):
            raise ValueError(f"The 'n_cores' parameter should be an integer"
                             f" but {self.n_cores} was given.")
        if self.select not in SELECT_LIST:
            raise ValueError(f"The 'select' parameter should be"
                             f" in {SELECT_LIST} but {self.select} was given.")
        # TODO self.start isinstance(TimeTree) and then this tree will be the starting tree
        #  Check for the mapping, and if not correct change accordingly
        if not (self.start in START_LIST or (isinstance(self.start, int) and self.start in range(len(trees)))
                or isinstance(self.start, TimeTreeSet)):
            raise ValueError(f"The 'start' parameter should be in {START_LIST} or "
                             f"an integer between 0 and {len(trees)-1} or"
                             f"a TimeTreeSet containing 1 tree!")

        # Computing/selecting the starting tree with the specified variation
        starting_tree = trees[0]  # Initializing with the case "first"

        if self.start == "last":
            starting_tree = trees[-1]
        elif self.start == "random":
            starting_tree = trees[random.sample(range(len(trees)), 1)[0]]
        elif isinstance(self.start, int):
            starting_tree = trees[self.start]
        elif self.start == "FM":
            starting_tree = frechet_mean(trees)
        elif isinstance(self.start, TimeTreeSet):
            if not len(self.start) == 1:
                raise ValueError("Starting with a TimeTreeSet only possible if it contains a single tree!")
                # maybe change this in the future
            if not len(trees[0]) == len(self.start[0]):
                raise ValueError(f"The given starting tree has different number of taxa! "
                                 f"({len(trees[0])} vs. {len(self.start[0])})")
            # todo check the mapping and change accordingly if nessecary
            starting_tree = self.start[0]  # Setting the start to the only tree in the set

        return getattr(_variations, self.variation)(trees=trees, n_cores=self.n_cores, select=self.select,
                                                    start=starting_tree, subsample_size=self.subsample_size)


if __name__ == '__main__':

    from timeit import default_timer as timer

    import numpy as np

    # d_name = "binary_K047_no_error_including_burnin"

    mystart = TimeTreeSet('/Users/larsberling/Desktop/CodingMA/Git/Summary/Simulations/31_800_005_0/conv_centroid.tree')

    myts = TimeTreeSet('/Users/larsberling/Desktop/CodingMA/Git/Summary/Simulations/31_800_005_0/31_800_005_0.trees')

    mycen = Centroid(start=mystart, variation="inc_sub")

    cen, sos = mycen.compute_centroid(myts)

    print(sos)

    # d_name = 'RSV2'
    #
    # for i in range(2):
    #     myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')
    #
    #     mycen = Centroid(start="FM", variation="inc_sub")
    #
    #     cen, sos = mycen.compute_centroid(myts)
    #
    #     print(sos)
