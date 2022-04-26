import arviz as az
import numpy as np
import pandas as pd
import xarray as xr

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

coda = importr("coda")
tracerer = importr("tracerer")


def coda_ess(data_list, chain_length, sampling_interval):
    coda_mcmc = coda.mcmc(FloatVector(list(data_list)), start=0, end=chain_length,
                          thin=sampling_interval)
    return coda.effectiveSize(coda_mcmc)[0]


def tracerer_ess(data_list, sampling_interval, **kwargs):
    return tracerer.calc_ess(FloatVector(list(data_list)), sample_interval=sampling_interval)[0]


def arviz_ess(data_list, **kwargs):
    return az.ess(np.asarray(data_list))


# todo this should take the distance measure as an argument, i.e. RF, PD or RNNI as options
#  First two from the RWTY original
# todo it should also be possible to change the return by giving median, mean or min as option
#  see JC paper for this
import random

def pseudo_ess(tree_set, chain_length, sampling_interval, **kwargs):
    # todo this is still WIP!
    ess = []
    # todo the range parameter should also be an option
    for _ in range(10):
        cur_focal_fix = random.randint(0, len(tree_set)-1)
        cur_distance_list = [t.fp_distance(tree_set[cur_focal_fix]) for t in tree_set]
        ess.append(coda_ess(cur_distance_list, chain_length=chain_length, sampling_interval=sampling_interval))
    return np.mean(ess)

def _ess_tracerer_rsample():
    # todo tracerer on a random sample of values
    tracerer = importr("tracerer")
    sample = np.random.randn(1000, 1)
    ess = tracerer.calc_ess(FloatVector(sample), sample_interval=1)

    print(ess)


def _multi_dim_ess():

    # todo this could be interesting for future applications on mds coordinates

    sample = np.random.randn(200, 3)
    for _ in range(2):
        sample = np.append(sample, sample, axis=0)
    sample = pd.DataFrame({"a":sample[:,0], "b": sample[:,1], "c": sample[:,2]})
    sample["chain"] = 0
    sample["draw"] = np.arange(len(sample), dtype=int)
    sample = sample.set_index(["chain", "draw"])
    xsample = xr.Dataset.from_dataframe(sample)

    dataset = az.InferenceData(posterior=xsample)
    print(az.ess(dataset))

    return 0
