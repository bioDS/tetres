import numpy as np
import pandas as pd


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
        ess.append(calc_ess(data_list=cur_distance_list))
    return np.median(ess)


def calc_ess(trace) -> float:
    """
    Returns the effective sample size (ESS) estimate of the input trace.

    :param trace: Trace of values, can be anything that is translatable via 'np.asarray()'
    :return: Effective sample size estimate of the input trace
    :rtype: float
    """
    # Reimplementation of ess as in the ASM package, translated from java
    trace = np.asarray(trace, dtype=np.float64)
    return len(trace) / (act(trace))


def act(trace, sample_interval: int = 1) -> np.float64:
    """
    Calcualtes the autocorrelation time for a given trace and a specific sampling interval

    :param trace: Trace of values, can be anything that is translatable via 'np.asarray()'
    :param sample_interval: integer specifying the interval between samples
    :return: autocorrelation time of independent samples
    :rtype: np.float64
    """
    MAX_LAG = 2000  # for efficiency
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
