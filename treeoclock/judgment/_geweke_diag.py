from treeoclock.trees.time_trees import TimeTreeSet

import numpy as np

import itertools


def geweke_diagnostic_summary():

    return 0


def geweke_diagnostic_distances(trees: TimeTreeSet, norm: bool=False):
    new_log_list = [1]
    distance_list = {f"{r},{s}": trees.fp_distance(r, s, norm=norm) ** 2 for r, s in
                     list(itertools.permutations(range(len(trees)), 2))}
    intersum_list = [1]
    list10 = [1]
    list40 = [1]

    for i in range(1, len(trees)):
        if i < 10:
            # Setting 10 to be the smallest tree set for which the value is actually computed
            new_log_list.append(1)
            intersum_list.append(1)
            list10.append(1)
            list40.append(1)
        else:
            sec10sum = 0
            last40sum = 0
            intersum = 0
            intersum_division = 0
            cur_10 = list(itertools.permutations(range(int(i * 0.1), int(i * 0.2)), 2))
            # cur_10 = list(itertools.permutations(range(int(i * 0.1)), 2))
            cur_40 = list(itertools.permutations(range(int(i * 0.6), i), 2))
            for r, s in cur_10:
                sec10sum += distance_list[f"{r},{s}"]
            for r, s in cur_40:
                last40sum += distance_list[f"{r},{s}"]
            for r in range(int(i * 0.1), int(i * 0.2)):
                # for r in range(int(i * 0.1)):
                for s in range(int(i * 0.6), i):
                    intersum_division += 1
                    intersum += distance_list[f"{r},{s}"]

            result = ((sec10sum / np.max([len(cur_10), 1])) - (intersum / intersum_division)) + (
                        (last40sum / len(cur_40)) - (intersum / intersum_division))
            new_log_list.append(result)
            intersum_list.append((sec10sum / np.max([len(cur_10), 1])) - (last40sum / len(cur_40)))
            list10.append(sec10sum / np.max([len(cur_10), 1]))
            list40.append(last40sum / len(cur_40))

    # todo this is still WIP
    return new_log_list, intersum_list, list10, list40


def geweke_diagnostic_variacne():

    return 0


