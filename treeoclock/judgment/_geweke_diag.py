from treeoclock.trees.time_trees import TimeTreeSet

import numpy as np

import itertools


def _check_percentage_input(first_range, last_percent):
    if len(first_range) > 2:
        raise Warning("More than two values given!")
    if first_range[0] > first_range[1]:
        raise ValueError("Given Range is Empty!")
    range_check = [True if 0 < x < 1 else False for x in first_range]
    if not range_check:
        raise ValueError("Given Range is not between 0 and 1!")
    if not 0 < last_percent < 1:
        raise ValueError("Values needs to be between 0 and 1!")
    if 1-last_percent < first_range[1]:
        raise Warning("The given ranges overlap!")


def geweke_diagnostic_distances(trees: TimeTreeSet, norm: bool = False, first_range=[0.1, 0.2], last_percent=0.4):

    _check_percentage_input(first_range, last_percent)

    new_log_list = [1]
    distance_list = {f"{r},{s}": trees.fp_distance(r, s, norm=norm) ** 2 for r, s in
                     list(itertools.permutations(range(len(trees)), 2))}

    for i in range(1, len(trees)):
        if i < 10:
            # Setting 10 to be the smallest tree set for which the value is actually computed
            new_log_list.append(1)
        else:
            # todo rename all the variables with more general names
            sec10sum = 0
            last40sum = 0
            intersum = 0
            intersum_division = 0
            cur_10 = list(itertools.permutations(range(int(i * first_range[0]), int(i * first_range[1])), 2))
            cur_40 = list(itertools.permutations(range(int(i * (1-last_percent)), i), 2))
            for r, s in cur_10:
                sec10sum += distance_list[f"{r},{s}"]
            for r, s in cur_40:
                last40sum += distance_list[f"{r},{s}"]
            for r in range(int(i * first_range[0]), int(i * first_range[1])):
                for s in range(int(i * (1-last_percent)), i):
                    intersum_division += 1
                    intersum += distance_list[f"{r},{s}"]

            result = np.absolute((sec10sum / np.max([len(cur_10), 1])) - (last40sum / len(cur_40)))
            result += np.absolute((sec10sum / np.max([len(cur_10), 1])) - (intersum / intersum_division))
            result += np.absolute((intersum / intersum_division) - (last40sum / len(cur_40)))
            new_log_list.append(np.sqrt(result))

    return new_log_list


def geweke_diagnostic_variance():
    return 0


def geweke_diagnostic_summary():
    return 0

