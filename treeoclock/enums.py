from enum import Enum
import numpy as np


def _median_ad(*args):
    return np.median(np.absolute(args - np.median(args)))
def _mean_ad(*args):
    return np.mean(np.absolute(args - np.mean(args)))


class Summary(Enum):
    FM = 'frechet_mean'
    CEN = 'centroid'
    TREE = 'tree'


class Avg_enum(Enum):
    MEAN = np.mean
    MEDIAN = np.median
    MEDIAN_AD = _median_ad
    MEAN_AD = _mean_ad
