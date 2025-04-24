import numpy as np
import pandas as pd

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
tracerer = importr("tracerer")


MAX_LAG = 2000  # for efficiency


def _ess_df(cmchain, chain_indeces, start=-1, end=-1):
    if start == -1:
        start = 0
    df = []
    end_internal = end
    for chain in chain_indeces:
        if end == -1:
            end_internal = len(cmchain[chain].trees) - 1
        for k in cmchain[chain].log_data.keys():
            if k != "Sample":
                df.append([k, cmchain[chain].get_ess(ess_key=k, lower_i=start, upper_i=end_internal), cmchain[chain].name])
        df.append(
            ["Pseudo_ESS_RNNI", cmchain[chain].get_pseudo_ess(sample_range=100, lower_i=start, upper_i=end_internal), cmchain[chain].name])
        df.append(
            ["Pseudo_ESS_RF", cmchain[chain].get_pseudo_ess(dist="rf", sample_range=100, lower_i=start, upper_i=end_internal), cmchain[chain].name])

    return pd.DataFrame(df, columns=["Key", "Value", "Chain"])


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
    max_lag = min(n, max_lag)
    # return the ESS or the number of samples if ESS overestimates the value
    return n / (-1 + (2 * np.sum(_autocorr_t(data_list, max_lag, trunc))))


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
    m = np.mean(data_list)
    # m = sum(data_list)/n
    # gamma_0 = sum([(y - m)**2 for y in data_list])
    gamma_0 = (n-1)*np.var(data_list)
    def gamma(k):
        nonlocal data_list
        return np.sum([(data_list[i] - m) * (data_list[i+k] - m) for i in range(0,n-k)])
    cor = list()
    for i in range(max_lag):
        c = gamma(i) / gamma_0
        if c < trunc:
            break
        cor.append(c)
    return cor


def pseudo_ess(tree_set, dist="rnni", sample_range=10, no_zero=False):
    # todo use the pairwise distance matrix if avaialble?!
    #  For this i will need to add upper and lower i boundaries for this function
    #  and then calculate with that!
    #  also the parameter of rf needs to be given to the pwd_matrix funciton!
    ess = []
    # samples = random.sample(range(len(tree_set)), sample_range)
    samples = np.linspace(0, len(tree_set)-1, num=sample_range, dtype=int)

    for s in samples:
        cur_focal_fix = tree_set[s]
        if dist == "rf":
            cur_distance_list = [int(cur_focal_fix.etree.robinson_foulds(t.etree)[0]) for t in tree_set]
        elif dist == "rnni":
            cur_distance_list = [t.fp_distance(cur_focal_fix) for t in tree_set]
        else:
            raise ValueError(f"Unkown distance given {dist}!")
        if no_zero:
            cur_distance_list = [d for d in cur_distance_list if d != 0]
        ess.append(autocorr_ess(data_list=cur_distance_list))
    return np.median(ess)


def calc_ess_beast2(trace, sample_interval=1):
    # Reimplementation of ess as in the ASM package, translated from java
    trace = np.asarray(trace, dtype=np.float64)
    return len(trace) / (act(trace, sample_interval) / sample_interval)


def act(trace, sample_interval=1):
    trace = np.asarray(trace, dtype=np.float64)
    n = len(trace)

    if n <= 1:
        return 0.0  # or np.nan depending on how you want to handle it

    square_lagged_sums = np.zeros(MAX_LAG)
    auto_correlation = np.zeros(MAX_LAG)
    total_sum = 0.0

    for i in range(n):
        total_sum += trace[i]
        mean = total_sum / (i + 1)

        sum1 = total_sum
        sum2 = total_sum

        max_lag_i = min(i + 1, MAX_LAG)
        for lag in range(max_lag_i):
            square_lagged_sums[lag] += trace[i - lag] * trace[i]
            auto_correlation[lag] = (
                square_lagged_sums[lag]
                - (sum1 + sum2) * mean
                + mean * mean * (i + 1 - lag)
            )
            auto_correlation[lag] /= (i + 1 - lag)

            sum1 -= trace[i - lag]
            sum2 -= trace[lag]

    max_lag = min(n, MAX_LAG)
    integral = 0.0

    for lag in range(max_lag):
        if lag == 0:
            integral = auto_correlation[0]
        elif lag % 2 == 0:
            if auto_correlation[lag - 1] + auto_correlation[lag] > 0:
                integral += 2.0 * (auto_correlation[lag - 1] + auto_correlation[lag])
            else:
                break

    return sample_interval * integral / auto_correlation[0] if auto_correlation[0] != 0 else np.inf


def tracer_ess_convenience(x):
    # Reimplementation of tracer ess, translated from cpp (convenience R package)
    x = np.asarray(x)
    samples = len(x)

    if samples <= 1:
        # edge cases would result in error.
        return 0.0

    max_lag = min(len(x)-1, MAX_LAG)

    mean = np.mean(x)
    gamma_stat = np.zeros(max_lag)

    for lag in range(max_lag):
        deltas = x[:samples - lag] - mean
        deltas_lagged = x[lag:] - mean
        gamma_stat[lag] = np.mean(deltas * deltas_lagged)

        if lag == 0:
            var_stat = gamma_stat[0]
        elif lag % 2 == 0:
            if gamma_stat[lag - 1] + gamma_stat[lag] > 0:
                var_stat += 2.0 * (gamma_stat[lag - 1] + gamma_stat[lag])
            else:
                max_lag = lag
                break
    act = var_stat / gamma_stat[0]
    ess = samples / act

    return ess
