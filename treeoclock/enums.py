from enum import Enum
import numpy as np
from collections import namedtuple


enum_entry = namedtuple("enum_entry", "name function")


def _median_ad(*args):
    return np.median(np.absolute(args - np.median(args)))
def _mean_ad(*args):
    return np.mean(np.absolute(args - np.mean(args)))


class Summary(Enum):
    FM = 'frechet_mean'
    CEN = 'centroid'
    TREE = 'tree'


class Avg_enum(Enum):
    MEAN = enum_entry("mean", np.mean)
    MEDIAN = enum_entry("median", np.median)
    MEDIAN_AD = enum_entry("median_ad", _median_ad)
    MEAN_AD = enum_entry("mean_ad", _mean_ad)
