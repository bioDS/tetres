from treeoclock.trees.time_trees import TimeTreeSet

import numpy as np
import itertools


def _check_percentage_input(first_range, last_percent):
    if len(first_range) > 2:
        raise ValueError("More than two values given!")
    if first_range[0] > first_range[1] or first_range[0] == first_range[1]:
        raise ValueError("Given Range is Empty!")
    range_check = [True if 0 < x < 1 else False for x in first_range]
    if not range_check:
        raise ValueError("Given Range is not between 0 and 1!")
    if not 0 < last_percent < 1:
        raise ValueError("Values needs to be between 0 and 1!")
    if 1 - last_percent < first_range[1]:
        raise ValueError("The given ranges overlap!")


def geweke_diagnostic_distances(trees: TimeTreeSet, norm: bool = False, first_range=[0.1, 0.2], last_percent=0.4):
    _check_percentage_input(first_range, last_percent)

    new_log_list = [1]  # Initialized list because the loop starts at 1

    distance_list = {f"{r},{s}": trees.fp_distance(r, s, norm=norm) ** 2
                 for r, s in list(itertools.permutations(range(len(trees)), 2))}

    for i in range(1, len(trees)):
        if i < 10:
            # Setting 10 to be the smallest tree set for which the value is actually computed
            new_log_list.append(1)
            # todo maybe instead check if either of the sets only contains a single tree and then append 1 or something
        else:
            first_set_sum = 0
            second_set_sum = 0
            intersum = 0
            intersum_division = 0
            first_set = list(itertools.permutations(range(int(i * first_range[0]), int(i * first_range[1])), 2))
            second_set = list(itertools.permutations(range(int(i * (1 - last_percent)), i), 2))
            for r, s in first_set:
                first_set_sum += distance_list[f"{r},{s}"]
            for r, s in second_set:
                second_set_sum += distance_list[f"{r},{s}"]
            for r in range(int(i * first_range[0]), int(i * first_range[1])):
                for s in range(int(i * (1 - last_percent)), i):
                    intersum_division += 1
                    intersum += distance_list[f"{r},{s}"]

            result = np.absolute((first_set_sum / np.max([len(first_set), 1])) - (second_set_sum / len(second_set)))
            result += np.absolute((first_set_sum / np.max([len(first_set), 1])) - (intersum / intersum_division))
            result += np.absolute((intersum / intersum_division) - (second_set_sum / len(second_set)))
            new_log_list.append(np.sqrt(result))

    return new_log_list


def geweke_diagnostic_variance():
    return 0


def geweke_diagnostic_summary():
    return 0
