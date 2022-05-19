from treeoclock.trees.time_trees import TimeTreeSet, findpath_distance

from treeoclock.judgment.posterior_analysis import calc_pw_distances, calc_pw_distances_two_sets

import random
import itertools
import numpy as np


def gelman_rubin_distance_diagnostic(tree_set1: TimeTreeSet, tree_set2: TimeTreeSet, norm: bool = False,
                                     samples: int = 100):
    # todo time the distance matrix computation and text remco how long it takes to do this
    distances_pwts1 = calc_pw_distances(tree_set1)
    distances_pwts2 = calc_pw_distances(tree_set2)
    distance_pwts1ts2 = calc_pw_distances_two_sets(tree_set1, tree_set2)

    # var_difference = []
    # in_sample_var_list = []
    # between_sample_var_list = []
    psrf_like = []

    for _ in range(samples):
        ts1_sample = random.randint(0, len(tree_set1) - 1)

        in_sample_var = np.sum(distances_pwts1[ts1_sample, :]) + np.sum(distances_pwts1[:, ts1_sample])

        between_sample_var = np.sum(distance_pwts1ts2[ts1_sample, :])

        # in_sample_var_list.append(in_sample_var / len(tree_set1))
        # between_sample_var_list.append(between_sample_var / len(tree_set2))
        # var_difference.append(abs((between_sample_var / len(tree_set2)) - (in_sample_var / len(tree_set1))))
        psrf_like.append(np.sqrt((between_sample_var / len(tree_set2)) / (in_sample_var / len(tree_set1))))

        # Same for a sample from the second tree set
        ts2_sample = random.randint(0, len(tree_set2) - 1)
        in_sample_var = np.sum(distances_pwts2[ts2_sample, :]) + np.sum(distances_pwts2[:, ts2_sample])
        between_sample_var = np.sum(distance_pwts1ts2[:, ts2_sample])
        psrf_like.append(np.sqrt((between_sample_var / len(tree_set1)) / (in_sample_var / len(tree_set2))))


    # return var_difference, in_sample_var_list, between_sample_var_list, psrf_like
    return psrf_like