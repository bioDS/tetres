from treeoclock.trees.time_trees import TimeTreeSet
from treeoclock.judgment._pairwise_distance_matrix import calc_pw_distances, calc_pw_distances_two_sets

import random
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def gelman_rubin_distance_diagnostic_plot(cMChain, samples: int = 100):

    figure, axis = plt.subplots(nrows=cMChain.m_MChains, ncols=cMChain.m_MChains, constrained_layout=True, figsize=[9, 7])

    # todo a cMChain initialize pwd_matrices using multiprocessing?

    # from distfit import distfit
    # dist = distfit()

    for i in range(cMChain.m_MChains-1):
        for j in range(i+1, cMChain.m_MChains):
            cur_psrf_like = gelman_rubin_distance_diagnostic_from_matrices(cMChain.pwd_matrix(i),
                                                                           cMChain.pwd_matrix(j),
                                                                           cMChain.pwd_matrix(i, j),
                                                                           samples=samples)

            sns.kdeplot(ax=axis[i, j], data=cur_psrf_like, x="PSRF_like", hue="Treeset", fill=True, legend=False,
                        common_norm=False)
            axis[i, j].axvline(x=np.mean(cur_psrf_like["PSRF_like"]), color="red")
            axis[i, j].axvline(x=np.mean(cur_psrf_like["PSRF_like"]) + np.std(cur_psrf_like["PSRF_like"]), color="red", linestyle="--")
            axis[i, j].axvline(x=np.mean(cur_psrf_like["PSRF_like"]) - np.std(cur_psrf_like["PSRF_like"]), color="red", linestyle="--")
            # print(np.var(cur_psrf_like["PSRF_like"]))
            # import scipy
            # x0, x1 = axis[i, j].get_xlim()
            # x_pdf = np.linspace(x0, x1, 100)
            # y_pdf = scipy.stats.norm.pdf(x_pdf, loc=1, scale=0.005)
            # sns.lineplot(x=x_pdf, y=y_pdf, color="red", ax=axis[i, j])

            # dist.fit_transform(cur_psrf_like["PSRF_like"])
            # # print(dist.summary)
            # dist.plot(title="", ax=axis[i, j])
            # axis[i, j].get_legend().remove()

            for label in axis[i, j].get_xticklabels():
                label.set_rotation(45)
            sns.boxplot(ax=axis[j, i], data=cur_psrf_like, x="Treeset", y="PSRF_like")

            # from statsmodels.stats.weightstats import ztest
            # zt = ztest(x1=cur_psrf_like[cur_psrf_like["Treeset"] == "TS1"]["PSRF_like"],
            #            x2=cur_psrf_like[cur_psrf_like["Treeset"] == "TS2"]["PSRF_like"],
            #            alternative="two-sided")
            # print(zt)
            # from scipy.stats import ks_2samp
            # ks = ks_2samp(data1=cur_psrf_like[cur_psrf_like["Treeset"] == "TS1"]["PSRF_like"],
            #          data2=cur_psrf_like[cur_psrf_like["Treeset"] == "TS2"]["PSRF_like"])
            # print(ks)

        cur_geweke = cMChain.MChain_list[i].compute_geweke_distances(index=i, name=cMChain.name, add=False)
        cur_sample = range(0, cMChain[i].chain_length + cMChain[i].tree_sampling_interval, cMChain[i].tree_sampling_interval)
        sns.lineplot(ax=axis[i, i], y=cur_geweke, x=cur_sample)
    # adding the last diagonal plot, because i will not run far enough
    cur_geweke = cMChain[-1].compute_geweke_distances(index=cMChain.m_MChains-1, name=cMChain.name, add=False)
    cur_sample = range(0, cMChain[-1].chain_length + cMChain[-1].tree_sampling_interval, cMChain[-1].tree_sampling_interval)
    sns.lineplot(ax=axis[cMChain.m_MChains-1, cMChain.m_MChains-1], y=cur_geweke, x=cur_sample)

        

    # Annotation of the plot
    pad = 5  # in points
    cols = [f"chain{m}" for m in range(cMChain.m_MChains)]
    # Setting the col names
    for ax, col in zip(axis[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords="axes fraction", textcoords="offset points",
                    size="large", ha="right", va="baseline")
    # Setting row names
    for ax, col in zip(axis[:, 0], cols):
        ax.annotate(col, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords="offset points",
                    size="large", ha="right", va="center")

    plt.savefig(fname=f"{cMChain.working_dir}/{cMChain.name}_grd_plot{'' if samples == 100 else f'_{samples}'}.png", format="png", bbox_inches="tight", dpi=800)


def gelman_rubin_distance_diagnostic_from_matrices(pwts1, pwts2, pwts1ts2,
                                     samples: int = 100):
    nt1 = pwts1.shape[0]
    nt2 = pwts2.shape[0]
    psrf_like = []

    if samples != "all":
        for _ in range(samples):
            ts1_sample = random.randint(0, nt1 - 1)
            in_sample_var = np.sum(pwts1[ts1_sample, :]) + np.sum(pwts1[:, ts1_sample])
            between_sample_var = np.sum(pwts1ts2[ts1_sample, :])
            psrf_like.append(["TS1", np.sqrt((between_sample_var / nt2) / (in_sample_var / nt1))])

            ts2_sample = random.randint(0, nt2 - 1)
            in_sample_var = np.sum(pwts2[ts2_sample, :]) + np.sum(pwts2[:, ts2_sample])
            between_sample_var = np.sum(pwts1ts2[:, ts2_sample])
            psrf_like.append(["TS2", np.sqrt((between_sample_var / nt1) / (in_sample_var / nt2))])

    elif samples == "all":
        for cur_sample in range(nt1):
            in_sample_var = np.sum(pwts1[cur_sample, :]) + np.sum(pwts1[:, cur_sample])
            between_sample_var = np.sum(pwts1ts2[cur_sample, :])
            psrf_like.append(["TS1", np.sqrt((between_sample_var / nt2) / (in_sample_var / nt1))])
        for cur_sample in range(nt2):
            in_sample_var = np.sum(pwts2[cur_sample, :]) + np.sum(pwts2[:, cur_sample])
            between_sample_var = np.sum(pwts1ts2[:, cur_sample])
            psrf_like.append(["TS2", np.sqrt((between_sample_var / nt1) / (in_sample_var / nt2))])
    else:
        raise ValueError(f"Unrecognizes samples {samples} given!")
    return pd.DataFrame(data=psrf_like, columns=["Treeset", "PSRF_like"])


def gelman_rubin_diagnostic_list(cMChain, i, j):
    # Will compute the gelman rubin like diagnostic for chains i and j, this will result in a 'trace'

    # todo check index i and j for validity

    # todo check that both chains have same sized tree set




    return 0
