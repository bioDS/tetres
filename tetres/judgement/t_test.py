import numpy as np
import scipy.stats

def estimated_sample_var(trees, center):
    # todo this can probably be sped up using multithreading/processing

    s = sum([center.fp_distance(t, norm=True)**2 for t in trees])
    return s / len(trees)


def t_test(trees1, trees2, centroid1, centroid2, _var_version="similar"):
    if _var_version not in ["similar", "unequal"]:
        raise ValueError(f"Unsupported _var_version value {_var_version}")
    # Calculate the two estimated pooled variances
    s1 = estimated_sample_var(trees1, centroid1)
    s2 = estimated_sample_var(trees2, centroid2)
    # Check that the two variances are similar enough i.e. 1/2 < s1/s2 < 2
    # accoding to t test with different variances (not Welch t-test)
    if _var_version == "unequal":
    # if not 1/2 < s1/s2 < 2:
        # raise ValueError(f"Variances are out of similar ratio {s1/s2}, test not appropriate!")
        s_unbias = np.sqrt((s1/len(trees1)) + (s2/len(trees2)))
        t_statistic = (centroid1.fp_distance(centroid2, norm=True)) / (s_unbias * np.sqrt((1 / len(trees1)) + (1 / len(trees2))))
        fs1 = s1/len(trees1)
        fs2 = s2/len(trees2)
        df = ((fs1 + fs2)**2)/((fs1**2/(len(trees1)-1)) + (fs2**2/(len(trees2)-1)))
        return t_statistic, scipy.stats.t.sf(abs(t_statistic),
                                             df=df)
    elif _var_version == "similar":
        s_pooled = np.sqrt((((len(trees1) -1)*s1)+((len(trees2) -1)*s2))/(len(trees1)+len(trees2)-2))
        t_statistic = (centroid1.fp_distance(centroid2, norm=True))/(s_pooled * np.sqrt((1/len(trees1))+(1/len(trees2))))
        # Returning the t_statistic and the p value
        return t_statistic, scipy.stats.t.sf(abs(t_statistic), df=(len(trees1)+len(trees2)-2))*2
