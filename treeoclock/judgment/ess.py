import arviz as az
import numpy as np
import pandas as pd
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt
import random

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


def pseudo_ess(ess_method, tree_set, chain_length, sampling_interval, dist="rnni", sample_range=10):
    ess = []
    samples = random.sample(range(len(tree_set)), sample_range)

    for cur_focal_fix_nbr in samples:
        cur_focal_fix = tree_set[cur_focal_fix_nbr]
        if dist == "rf":
            cur_distance_list = [cur_focal_fix.etree.robinson_foulds(t.etree)[0] for t in tree_set]
        elif dist == "rnni":
            cur_distance_list = [t.fp_distance(cur_focal_fix) for t in tree_set]
        else:
            raise ValueError(f"Unkown distance given {dist}!")
        ess.append(globals()[f"{ess_method}_ess"](data_list=cur_distance_list, chain_length=chain_length, sampling_interval=sampling_interval))
    # todo this can be changed to median or mean, look at JC paper min is best and median also good
    #  however mean seems to be not a good choice,
    #  in the test case min is the only one not overestimating the ESS
    return np.median(ess)  # median seems to be better!


def ess_stripplot(cmchain, ess_method):
    burnin = 0
    df = []
    for chain in cmchain:
        for k in chain.log_data.keys():
            if k != "Sample":
                df.append([k, chain.get_ess(ess_key=k, ess_method=ess_method), chain.name])
        df.append(["Pseudo_ESS_RNNI", chain.get_pseudo_ess(ess_method=ess_method, sample_range=100), chain.name])
        df.append(["Pseudo_ESS_RF", chain.get_pseudo_ess(ess_method=ess_method, dist="rf", sample_range=100), chain.name])

    df = pd.DataFrame(df, columns=["Key", "Value", "Chain"])
    ax = sns.stripplot(data=df, x="Key", y="Value", hue="Chain")

    for label in ax.get_xticklabels():
        label.set_rotation(90)
    plt.xlabel("")
    plt.ylabel("ESS")
    plt.suptitle(f"Number of samples: {len(cmchain[0].trees)}, burn-in={burnin*100}%")
    plt.tight_layout()
    # plt.show()

    plt.savefig(fname=f"{cmchain.working_dir}/plots/ess_{ess_method}_comparison.png", dpi=400)
    plt.clf()
    plt.close()


# todo implement a pseudo ess value using the centroid or fm tree instead of a random fixed focal tree


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
