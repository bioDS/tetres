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


from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary.frechet_mean import frechet_mean


# todo missing default values!
def geweke_diagnostic_focal_tree(trees: TimeTreeSet, focal_tree, norm, kind, first_range, last_percent):
    _check_percentage_input(first_range, last_percent)

    new_log_list = [1]
    for i in range(1, len(trees)):
        if focal_tree is "FM":
            if i < 10:  # todo temporary
                # Setting 10 to be the smallest tree set for which the value is actually computed
                new_log_list.append(1)
            else:
                # todo rename the variables, see above naming
                sec10 = trees[int(i * 0.1):int(i * 0.2)]
                last40 = trees[int(i * 0.6):i]

                if kind == "default":
                    var10 = compute_sos_mt(frechet_mean(sec10), sec10, norm=norm) / len(sec10)
                    var40 = compute_sos_mt(frechet_mean(last40), last40, norm=norm) / len(last40)
                    new_log_list.append(abs(var10 - var40))
                elif kind == "crossed":
                    var10_in40 = compute_sos_mt(frechet_mean(sec10), last40, norm=norm) / len(last40)
                    var40_in10 = compute_sos_mt(frechet_mean(last40), sec10, norm=norm) / len(sec10)
                    new_log_list.append(abs(var10_in40 - var40_in10))
                elif kind == "doublecrossed":
                    # compares the variation of two different trees for the same set
                    fm10 = frechet_mean(sec10)
                    fm40 = frechet_mean(last40)
                    var10 = compute_sos_mt(fm10, sec10, norm=norm) / len(sec10)
                    var40 = compute_sos_mt(fm40, last40, norm=norm) / len(last40)
                    var10_in40 = compute_sos_mt(fm10, last40, norm=norm) / len(last40)
                    var40_in10 = compute_sos_mt(fm40, sec10, norm=norm) / len(sec10)
                    new_log_list.append(abs(var10_in40 - var10) + abs(var40_in10 - var40))
                elif kind == "crosscompare":
                    # compares the variation of two trees in one set
                    fm10 = frechet_mean(sec10)
                    fm40 = frechet_mean(last40)
                    var10 = compute_sos_mt(fm10, sec10, norm=norm) / len(sec10)
                    var40 = compute_sos_mt(fm40, last40, norm=norm) / len(last40)
                    var10_in40 = compute_sos_mt(fm10, last40, norm=norm) / len(last40)
                    var40_in10 = compute_sos_mt(fm40, sec10, norm=norm) / len(sec10)
                    new_log_list.append(abs(var40_in10 - var10) + abs(var10_in40 - var40))
                else:
                    raise ValueError(f"The given kind {kind} is not recognized!")

        elif focal_tree is "Centroid":
            raise ValueError("Not yet implemented!")
        elif focal_tree is "random":
            raise ValueError("Not yet implemented!")
            # todo this actually does not need to have a summary tree? just use a randomly fixed tree or two randomly fixed trees?
            #  compare the results
    return new_log_list

