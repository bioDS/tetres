from treeoclock.trees.time_trees import TimeTreeSet, findpath_distance

from treeoclock.judgment.posterior_analysis import calc_pw_distances, calc_pw_distances_two_sets

import random
import itertools

def gelman_rubin_distance_diagnostic(tree_set1: TimeTreeSet, tree_set2: TimeTreeSet, norm: bool = False, samples: int = 100):

    print("Starting1")
    distances_pwts1 = calc_pw_distances(tree_set1)
    # distances_pwts1 = {f"{r},{s}": tree_set1.fp_distance(r, s, norm=norm) ** 2
    #                  for r, s in list(itertools.permutations(range(len(tree_set1)), 2))}
    # distances_pwts2 = {f"{r},{s}": tree_set2.fp_distance(r, s, norm=norm) ** 2
    #                  for r, s in list(itertools.permutations(range(len(tree_set2)), 2))}
    print("Starting12")
    # todo this should also be possible with a special mutliprocessing c function that takes two treesets as input!
    distance_pwts1ts2 = calc_pw_distances_two_sets(tree_set1, tree_set2)
    # distance_pwt1ts2 = {f"{r},{s}": findpath_distance(tree_set1[r].ctree, tree_set2[s].ctree, norm=norm) ** 2
    #                  for r in range(len(tree_set1)) for s in range(len(tree_set2))}

    var_difference = []

    for _ in range(samples):
        ts1_sample = random.randint(0, len(tree_set1) - 1)

        in_sample_var = 0
        between_sample_var = 0
        for i in range(len(tree_set1)):
            if i > ts1_sample:
                # in_sample_var += distances_pwts1[f"{ts1_sample},{i}"]
                in_sample_var += distances_pwts1[ts1_sample, i]
            elif i < ts1_sample:
                # in_sample_var += distances_pwts1[f"{i},{ts1_sample}"]
                in_sample_var += distances_pwts1[ts1_sample, i]

            # between_sample_var += distance_pwt1ts2[f"{ts1_sample},{i}"]
            between_sample_var += distance_pwts1ts2[ts1_sample, i]
        var_difference.append(abs((between_sample_var/len(tree_set2)) - (in_sample_var/len(tree_set1))))

        # ts2_sample = tree_set2[random.randint(0, len(tree_set2) - 1)]
    return var_difference
