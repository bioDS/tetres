from treeoclock.summary import _variations
from treeoclock.summary._constants import VAR_LIST, SELECT_LIST, START_LIST
from treeoclock.trees.time_trees import TimeTreeSet


# TODO Tests and Documentation


class Centroid:
    """
    Class to compute the centroid algorithm with a set of parameters
    """

    # TODO Default start is going to be FM, not implemented yet
    def __init__(self, variation="greedy", n_cores=None, select='random', start='last'):
        self.variation = variation  # Which centroid variation to compute
        self.n_cores = n_cores  # How many cores to use whenever Multiprocessing is/will be used
        self.select = select  # Specifying which tree to choose in case of multiple options with the same quality
        self.start = start  # Specifying the starting position for the centroid algorithm

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
        if not (self.start in START_LIST or (isinstance(self.start, int) and self.start in range(len(trees)))):
            raise ValueError(f"The 'start' parameter should be in {START_LIST} or "
                             f"an integer between 0 and {len(trees)-1}!")

        # TODO add a starting parameter to the class, FM, eFM, index, ...

        # Curently starting computation at the last tree of the set
        return getattr(_variations, self.variation)(trees=trees, n_cores=self.n_cores, select=self.select, start=trees[-1])


if __name__ == '__main__':
    d_name = "40_800_005"

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    mycen = Centroid()

    from timeit import default_timer as timer

    import numpy as np

    times = []
    for _ in range(1):
        s = timer()
        cen = mycen.compute_centroid(myts)
        times.append(timer()-s)
    print(np.mean(times))

    from treeoclock.summary.compute_sos import compute_sos_mt

    print(compute_sos_mt(t=myts[-1], trees=myts))
    print(compute_sos_mt(t=cen, trees=myts))