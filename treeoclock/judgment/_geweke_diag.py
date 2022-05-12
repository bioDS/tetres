import numpy as np
import itertools
import random

from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary.frechet_mean import frechet_mean
from treeoclock.trees.time_trees import TimeTreeSet
from treeoclock.summary.centroid import Centroid


def _treeset_centroid(trees: TimeTreeSet):
    # Currently the only way to change the centroid variation for geweke diagnostic is to change this piece of code
    return Centroid().compute_centroid(trees)[0]


def _treeset_random(trees: TimeTreeSet):
    return trees[random.randint(0, len(trees) - 1)]


_focal_tree_functions = {"FM": frechet_mean, "Centroid": _treeset_centroid, "Random": _treeset_random}
_geweke_kind = {"default", "crossed", "doublecrossed", "crosscompare"}


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

    new_log_list = []  # Initialized list because the loop starts at 1

    distance_list = {f"{r},{s}": trees.fp_distance(r, s, norm=norm) ** 2
                     for r, s in list(itertools.permutations(range(len(trees)), 2))}

    for i in range(0, len(trees)):
        if i < 10:
            # Setting 10 to be the smallest tree set for which the value is actually computed
            new_log_list.append(1)
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


def geweke_diagnostic_focal_tree(trees: TimeTreeSet, focal_tree: str = "FM", norm: bool = False, kind: str = "default",
                                 first_range=[0.1, 0.2], last_percent=0.4):
    _check_percentage_input(first_range, last_percent)

    if focal_tree not in _focal_tree_functions:
        raise ValueError("Given Focal Tree not accepted!")
    if kind not in _geweke_kind:
        raise ValueError("Given kind not accepted!")

    new_log_list = []
    for i in range(0, len(trees)):
        # if focal_tree is "FM":
        if i < 10:
            # Setting 10 to be the smallest tree set for which the value is actually computed
            new_log_list.append(1)
        else:
            first_set = trees[int(i * 0.1):int(i * 0.2)]
            second_set = trees[int(i * 0.6):i]

            first_focal = _focal_tree_functions[focal_tree](first_set)
            second_focal = _focal_tree_functions[focal_tree](second_set)

            if kind == "default":
                variance_first_set = compute_sos_mt(first_focal, first_set, norm=norm) / len(first_set)
                variance_second_set = compute_sos_mt(second_focal, second_set, norm=norm) / len(second_set)
                new_log_list.append(abs(variance_first_set - variance_second_set))
            elif kind == "crossed":
                variance_firstfocal_second_set = compute_sos_mt(first_focal, second_set, norm=norm) / len(second_set)
                variance_secondfocal_first_set = compute_sos_mt(second_focal, first_set, norm=norm) / len(first_set)
                new_log_list.append(abs(variance_firstfocal_second_set - variance_secondfocal_first_set))
            elif kind == "doublecrossed":
                # compares the variation of two different trees for the same set
                variance_first_set = compute_sos_mt(first_focal, first_set, norm=norm) / len(first_set)
                variance_second_set = compute_sos_mt(second_focal, second_set, norm=norm) / len(second_set)
                variance_firstfocal_second_set = compute_sos_mt(first_focal, second_set, norm=norm) / len(second_set)
                variance_secondfocal_first_set = compute_sos_mt(second_focal, first_set, norm=norm) / len(first_set)
                new_log_list.append(abs(variance_firstfocal_second_set - variance_first_set) +
                                    abs(variance_secondfocal_first_set - variance_second_set))
            elif kind == "crosscompare":
                # compares the variation of two trees in one set
                variance_first_set = compute_sos_mt(first_focal, first_set, norm=norm) / len(first_set)
                variance_second_set = compute_sos_mt(second_focal, second_set, norm=norm) / len(second_set)
                variance_firstfocal_second_set = compute_sos_mt(first_focal, second_set, norm=norm) / len(second_set)
                variance_secondfocal_first_set = compute_sos_mt(second_focal, first_set, norm=norm) / len(first_set)
                new_log_list.append(abs(variance_secondfocal_first_set - variance_first_set) +
                                    abs(variance_firstfocal_second_set - variance_second_set))
            else:
                raise ValueError(f"The given kind {kind} is not recognized!")
    return new_log_list
