from treeoclock.trees.time_trees import TimeTreeSet
from treeoclock.judgment._pairwise_distance_matrix import calc_pw_distances, calc_pw_distances_two_sets

import random
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def gelman_rubin_distance_diagnostic_plot(cMChain, samples: int = 100):
    figure, axis = plt.subplots(nrows=cMChain.m_MChains, ncols=cMChain.m_MChains, constrained_layout=True,
                                figsize=[9, 7])

    # todo a cMChain initialize pwd_matrices using multiprocessing?

    # from distfit import distfit
    # dist = distfit()

    for i in range(cMChain.m_MChains - 1):
        for j in range(i + 1, cMChain.m_MChains):
            cur_psrf_like = gelman_rubin_distance_diagnostic_from_matrices(cMChain.pwd_matrix(i),
                                                                           cMChain.pwd_matrix(j),
                                                                           cMChain.pwd_matrix(i, j),
                                                                           samples=samples)

            sns.kdeplot(ax=axis[i, j], data=cur_psrf_like, x="PSRF_like", hue="Treeset", fill=True, legend=False,
                        common_norm=False)
            axis[i, j].axvline(x=np.mean(cur_psrf_like["PSRF_like"]), color="red")
            axis[i, j].axvline(x=np.mean(cur_psrf_like["PSRF_like"]) + np.std(cur_psrf_like["PSRF_like"]), color="red",
                               linestyle="--")
            axis[i, j].axvline(x=np.mean(cur_psrf_like["PSRF_like"]) - np.std(cur_psrf_like["PSRF_like"]), color="red",
                               linestyle="--")
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
        cur_sample = range(0, cMChain[i].chain_length + cMChain[i].tree_sampling_interval,
                           cMChain[i].tree_sampling_interval)
        sns.lineplot(ax=axis[i, i], y=cur_geweke, x=cur_sample)
    # adding the last diagonal plot, because i will not run far enough
    cur_geweke = cMChain[-1].compute_geweke_distances(index=cMChain.m_MChains - 1, name=cMChain.name, add=False)
    cur_sample = range(0, cMChain[-1].chain_length + cMChain[-1].tree_sampling_interval,
                       cMChain[-1].tree_sampling_interval)
    sns.lineplot(ax=axis[cMChain.m_MChains - 1, cMChain.m_MChains - 1], y=cur_geweke, x=cur_sample)

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

    plt.savefig(fname=f"{cMChain.working_dir}/{cMChain.name}_grd_plot{'' if samples == 100 else f'_{samples}'}.png",
                format="png", bbox_inches="tight", dpi=800)
    plt.clf()


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


def gelman_rubin_trace_with_cutoff(cmchain, i, j, sample_from, threshold_percentage):
    dmi = cmchain.pwd_matrix(i)
    dmj = cmchain.pwd_matrix(j)
    dmij = cmchain.pwd_matrix(i, j)

    if dmi.shape != dmj.shape:
        raise ValueError("Treesets have different sizes!")

    df = []

    cutoff = {k: -1 for k in threshold_percentage}
    consecutive = 0

    for cur_sample in range(1, dmi.shape[0]):
        # starting at 10 to use add and because small treesets don't make so much sense
        i_in_var = np.sum(dmi[cur_sample, :cur_sample]) + np.sum(dmi[:cur_sample, cur_sample])
        i_between_var = np.sum(dmij[cur_sample, :cur_sample])
        psrf_like = np.sqrt((i_between_var / cur_sample) / (i_in_var / cur_sample))
        psrf_like = [psrf_like]
        for x in range(int(cur_sample * (1 - sample_from)), cur_sample - 1):
            x_in_var = np.sum(dmi[x, :cur_sample]) + np.sum(dmi[:cur_sample, x])
            x_between_var = np.sum(dmij[x, :cur_sample])
            psrf_like.append(np.sqrt((x_between_var / cur_sample) / (x_in_var / cur_sample)))
        psrf_like = np.median(psrf_like)
        df.append([cur_sample, psrf_like, f"chain_{i}"])

        j_in_var = np.sum(dmj[cur_sample, :cur_sample]) + np.sum(dmj[:cur_sample, cur_sample])
        j_between_var = np.sum(dmij[:cur_sample, cur_sample])
        psrf_like = np.sqrt((j_between_var / cur_sample) / (j_in_var / cur_sample))
        psrf_like = [psrf_like]
        for x in range(int(cur_sample * (1 - sample_from)), cur_sample - 1):
            x_in_var = np.sum(dmj[x, :cur_sample]) + np.sum(dmj[:cur_sample, x])
            x_between_var = np.sum(dmij[:cur_sample, x])
            psrf_like.append(np.sqrt((x_between_var / cur_sample) / (x_in_var / cur_sample)))
        psrf_like = np.median(psrf_like)
        df.append([cur_sample, psrf_like, f"chain_{j}"])

        if 0.99 < df[-1][1] < 1.01 and 0.99 < df[-2][1] < 1.01:
            consecutive += 1
            for k in cutoff.keys():
                if cutoff[k] == -1 and consecutive >= int(k*cur_sample):
                    cutoff[k] = cur_sample
        else:
            consecutive = 0

    return pd.DataFrame(df, columns=["Sample", "PSRF", "Chain"]), cutoff


def gelman_rubin_trace_plot(cmchain, i, j):
    # Will compute the gelman rubin like diagnostic for chains i and j, this will result in a 'trace'

    sample_from = [0, 0.1, 0.5, 0.75, 1, 0.25]
    sample_from = np.sort(sample_from)
    threshold_percentage = [0, 0.1, 0.25, 0.5, 0.75, 1.0]
    threshold_percentage = np.sort(threshold_percentage)

    figure, axis = plt.subplots(ncols=len(sample_from), nrows=len(threshold_percentage),
                                constrained_layout=True,
                                # figsize=[9, 7],
                                squeeze=False, sharex=True, sharey=True)

    for col in range(len(sample_from)):
        df, cutoff = gelman_rubin_trace_with_cutoff(cmchain, i, j, sample_from=sample_from[col],
                                                    threshold_percentage=threshold_percentage)
        for row in range(len(threshold_percentage)):
            sns.lineplot(data=df, x="Sample", y="PSRF", hue="Chain", alpha=0.5, ax=axis[row, col], legend=False)
            # axis[row, col].axhline(y=1.01, linestyle="--", color="red")
            # axis[row, col].axhline(y=0.99, linestyle="--", color="red")

            axis[row, col].set_xticks = set(df["Sample"])
            # todo the ytick could maybe be adapted
            #  for cases where the value is above, at 1.2 or 1.5 for diverging datasets?
            # axis[row, col].set_ylim([0.8, 1.2])
            # axis[row, col].set_yticks = np.arange(0.8, 1.2, 0.05)

            axis[row, col].axvline(x=cutoff[threshold_percentage[row]], color="red")
            axis[row, col].text(x=cutoff[threshold_percentage[row]] + 0.1, y=1.1,
                                s=f'{cutoff[threshold_percentage[row]]}',
                                fontsize=8, zorder=20)

            axis[row, col].set_ylabel(f"{threshold_percentage[row]}", color="green")
            axis[row, col].set_xlabel(f"{sample_from[col]}", color="blue")
    figure.supylabel("Threshold Time (%)", color="green")
    figure.supxlabel("Sample from last x% of trees", color="blue")

    # plt.axhline(y=1.05, linestyle="--", color="orange")
    # plt.axhline(y=0.95, linestyle="--", color="orange")
    # 
    # plt.axhline(y=1.1, linestyle="--", color="yellow")
    # plt.axhline(y=0.9, linestyle="--", color="yellow")

    plt.savefig(fname=f"{cmchain.working_dir}/{cmchain.name}_{i}-{j}_grd_singlevalue_evaluation.png",
                format="png", bbox_inches="tight", dpi=800)
    plt.clf()

