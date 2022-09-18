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
    # todo use the pairwise distance matrix if avaialble?! For this i will need to add upper and lower i boundaries for this function and then calculate with that!
    #  also the parameter of rf needs to be given to the pwd_matrix funciton!
    ess = []
    # samples = random.sample(range(len(tree_set)), sample_range)
    samples = np.linspace(0, len(tree_set)-1, num=sample_range, dtype=int)

    for s in samples:
        cur_focal_fix = tree_set[s]
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

    df = _ess_df(cmchain, chain_indeces=range(cmchain.m_MChains), ess_method=ess_method)

    ax = sns.stripplot(data=df, x="Key", y="Value", hue="Chain")

    for label in ax.get_xticklabels():
        label.set_rotation(90)
    plt.xlabel("")
    plt.ylabel("ESS")
    plt.suptitle(f"Number of samples: {len(cmchain[-1].trees)}, burn-in={burnin*100}%")
    plt.tight_layout()
    # plt.show()

    plt.savefig(fname=f"{cmchain.working_dir}/plots/{cmchain.name}_{ess_method}_ess_comparison.png", dpi=400)
    plt.clf()
    plt.close()


def _ess_df(cmchain, chain_indeces, ess_method, start=-1, end=-1):
    if start == -1:
        start = 0
    df = []
    end_internal = end
    for chain in chain_indeces:
        if end == -1:
            end_internal = len(cmchain[chain].trees) - 1
        for k in cmchain[chain].log_data.keys():
            if k != "Sample":
                df.append([k, cmchain[chain].get_ess(ess_key=k, ess_method=ess_method, lower_i=start, upper_i=end_internal), cmchain[chain].name])
        df.append(
            ["Pseudo_ESS_RNNI", cmchain[chain].get_pseudo_ess(ess_method=ess_method, sample_range=100, lower_i=start, upper_i=end_internal), cmchain[chain].name])
        df.append(
            ["Pseudo_ESS_RF", cmchain[chain].get_pseudo_ess(ess_method=ess_method, dist="rf", sample_range=100, lower_i=start, upper_i=end_internal), cmchain[chain].name])

    return pd.DataFrame(df, columns=["Key", "Value", "Chain"])


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


def autocorr_ess(data_list, max_lag=2000, trunc=0.05):
    """
    Effective Sample Size calculated using autocorrelation.
    
    Effective Sample Size (ESS) calculated using the autocorrelation method.
    The ESS is defined as: `n / tau`, where `tau = 1 + 2 sum_{k=1}^{\infty} tho_k(theta)`.
    An estimate of `tau` is used where the infinite sum is truncated at lag `k` where:
    `rho_k < 0.05` (see e.g., van Dyk and Park (2011) or Kass et al. (1988))
    
    Parameters:
        x (list): a list of numeric values
        max_lag (int): maximum calculated lag of the autocorrelation function
        trunc (double): stop calculating the autocorrelation at this value
    Returns:
        ess (double): estimate of the Effective Sample Size
    """
    n = len(data_list)
    max_lag = min(n - 1, max_lag)
    return n / (-1 + 2 * sum(_autocorr_t(data_list, max_lag, trunc)))


def _autocorr_t(data_list, max_lag=2000, trunc=0.05):
    """
    Calculate lag-k autocorrelation.
    
    Calculate the lag-k autocorrelation ussing the gamma function (see Duerre et al. 2014)
    
    The correlation coefficient `rho_k` for the lag `k` is calculated using the gamma function:
    `tho_k = gamma(k) / gamma(0)`, where `gamma(k) = sum_{j=1}^{n-k} (x_j - m) (x_{j+1}) - m)`.
    Note that instead of the mean for the partial series, a sample mean is used.
    
    Parameters:
        x (list): a list of numeric values
        max_lag (int): maximum calculated lag of the autocorrelation function
        trunc (double): stop calculating the atucorrelation function at this value
    Returns:
        cor (list): truncated autocorrelation series
    """
    n = len(data_list)
    m = sum(data_list)/n
    gamma_0 = sum([(y - m)**2 for y in data_list])
    def gamma(k):
        return sum([ (data_list[i] - m) * (data_list[i+k] - m) for i in range(0,n-k) ])
    cor = list()
    for i in range(max_lag):
        c = gamma(i) / gamma_0
        cor.append(c)
        if(c < trunc):
            break
    return cor
