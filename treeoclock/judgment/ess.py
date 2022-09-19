import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def ess_stripplot(cmchain):
    # todo why/where do i need this?
    burnin = 0

    df = _ess_df(cmchain, chain_indeces=range(cmchain.m_MChains))

    ax = sns.stripplot(data=df, x="Key", y="Value", hue="Chain")

    for label in ax.get_xticklabels():
        label.set_rotation(90)
    plt.xlabel("")
    plt.ylabel("ESS")
    plt.suptitle(f"Number of samples: {len(cmchain[-1].trees)}, burn-in={burnin*100}%")
    plt.tight_layout()
    # plt.show()

    plt.savefig(fname=f"{cmchain.working_dir}/plots/{cmchain.name}_ess_comparison.png", dpi=400)
    plt.clf()
    plt.close()


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
    max_lag = min(n - 1, max_lag)
    # return the ESS or the number of samples if ESS overestimates the value
    return min(n / (-1 + 2 * sum(_autocorr_t(data_list, max_lag, trunc))), n)


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


def pseudo_ess(tree_set, dist="rnni", sample_range=10):
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
        ess.append(autocorr_ess(data_list=cur_distance_list))
    return np.median(ess)
