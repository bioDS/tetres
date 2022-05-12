from treeoclock.summary import _variations
from treeoclock.summary._constants import SELECT_LIST, START_LIST
from treeoclock.trees.time_trees import TimeTreeSet
from treeoclock.summary.frechet_mean import frechet_mean, frechet_mean_sort

import os
import random
import warnings

# TODO Tests and Documentation


class Centroid:
    """
    Class to compute the centroid algorithm with a set of parameters
    """

    def __init__(self, variation="greedy_omp", n_cores=None, select='random', start='FM', subsample_size=200,
                 tree_log_file="", max_iterations=None):
        self.variation = variation  # Which centroid variation to compute
        self.n_cores = n_cores  # How many cores to use whenever Multiprocessing is/will be used
        self.select = select  # Specifying which tree to choose in case of multiple options with the same quality
        self.start = start  # Specifying the starting position for the centroid algorithm
        self.subsample_size = subsample_size  # Size of the subsamples for some of the centroid variations
        self.tree_log_file = tree_log_file  # If a file is given, the trace of trees will be output to this file
        self.max_iterations = max_iterations

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
        elif self.start == "sortFM":
            starting_tree = frechet_mean_sort(trees)
        elif isinstance(self.start, TimeTreeSet):
            if not len(trees[0]) == len(self.start[0]):
                raise ValueError(f"The given starting tree has different number of taxa! "
                                 f"({len(trees[0])} vs. {len(self.start[0])})")
            if not trees.map == self.start.map:
                self.start.change_mapping(new_map=trees.map)

            starting_tree = self.start[0]  # Setting the start to the first tree in the given set
            # todo maybe change the index with parameter at some point in the future

        if self.tree_log_file:
            # will write a log_file
            if not os.path.exists(os.path.dirname(self.tree_log_file)):
                raise FileNotFoundError('Given Path to logfile does not exist!')
            if os.path.exists(self.tree_log_file):
                warnings.warn(UserWarning('Given Tree Log File exists already and will be overwritten!'))
            with open(self.tree_log_file, "w+") as log:
                log.write(f"#NEXUS\n\nBegin taxa;\n\tDimensions ntax={len(trees[0])};\n\t\tTaxlabels\n")
                for k in sorted(trees.map.keys()):
                    log.write(f"\t\t\t{trees.map[k]}\n")
                log.write("\t\t\t;\nEnd;\n")

                log.write("Begin trees;\n\tTranslate\n")
                for k in sorted(trees.map.keys()):
                    log.write(f"\t\t{k} {trees.map[k]}")
                    if k is not len(trees.map.keys()):
                        log.write(',')
                    log.write('\n')
                log.write(";\n")

        return getattr(_variations, self.variation)(trees=trees, n_cores=self.n_cores, select=self.select,
                                                    start=starting_tree, subsample_size=self.subsample_size,
                                                    tree_log_file=self.tree_log_file, max_iterations=self.max_iterations)
