import gc
import os
import random
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from treeoclock.judgment.ess import _ess_df


def gelman_rubin_distance_diagnostic_plot(cMChain, samples: int = 100):
    figure, axis = plt.subplots(nrows=cMChain.m_MChains, ncols=cMChain.m_MChains, constrained_layout=True,
                                figsize=[9, 7])

    # initializing the scaled axis for the upper part of the plot:
    scaled_x_axis = np.arange(0.8, 1.7, 0.1)
    scaled_y_axis = np.arange(0, 120, 20)

    for i in range(cMChain.m_MChains - 1):
        for j in range(i + 1, cMChain.m_MChains):
            cur_psrf_like = gelman_rubin_distance_diagnostic_from_matrices(cMChain.pwd_matrix(i),
                                                                           cMChain.pwd_matrix(j),
                                                                           cMChain.pwd_matrix(i, j),
                                                                           samples=samples)

            # upper triangle of plot
            sns.kdeplot(ax=axis[i, j], data=cur_psrf_like, x="PSRF_like", hue="Treeset", fill=True, legend=False,
                        common_norm=False)
            axis[i, j].set_xticks(scaled_x_axis)
            axis[i, j].set_yticks(scaled_y_axis)
            axis[i, j].set_xlim(np.min(scaled_x_axis), np.max(scaled_x_axis))
            axis[i, j].set_ylim(np.min(scaled_y_axis), np.max(scaled_y_axis))
            for label in axis[i, j].get_xticklabels():
                label.set_rotation(90)

            # lower triangle of plot
            sns.kdeplot(ax=axis[j, i], data=cur_psrf_like, x="PSRF_like", hue="Treeset", fill=True, legend=False,
                        common_norm=False)
            # todo maybe add a line at 1 (if in the xtick range) indicating ideal and then the actual mean
            axis[j, i].axvline(x=np.mean(cur_psrf_like["PSRF_like"]), color="red")
            axis[j, i].axvline(x=np.mean(cur_psrf_like["PSRF_like"]) + np.std(cur_psrf_like["PSRF_like"]), color="red",
                               linestyle="--")
            axis[j, i].axvline(x=np.mean(cur_psrf_like["PSRF_like"]) - np.std(cur_psrf_like["PSRF_like"]), color="red",
                               linestyle="--")

            for label in axis[j, i].get_xticklabels():
                label.set_rotation(45)
            axis[j, i].set_xticks(np.linspace(axis[j, i].get_xticks().min(),
                                              axis[j, i].get_xticks().max(),
                                              5))

            # sns.boxplot(ax=axis[j, i], data=cur_psrf_like, x="Treeset", y="PSRF_like")

        # putting geweke on diagonal of the plot
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

    # plt.show()
    plt.savefig(
        fname=f"{cMChain.working_dir}/plots/{cMChain.name}_grd_plot{'' if samples == 100 else f'_{samples}'}.png",
        format="png", bbox_inches="tight", dpi=800)
    plt.clf()
    plt.close("all")
    gc.collect()


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


def gelman_rubin_trace_with_cutoff(cmchain, i, j, sample_from, threshold_percentage, ess=0, pess_range=100):
    dmi = cmchain.pwd_matrix(i)
    dmj = cmchain.pwd_matrix(j)
    dmij = cmchain.pwd_matrix(i, j)

    if dmi.shape != dmj.shape:
        raise ValueError("Treesets have different sizes!")

    df = []

    cutoff_start = {k: -1 for k in threshold_percentage}
    cutoff_end = {k: -1 for k in threshold_percentage}
    consecutive = 0

    for cur_sample in range(dmi.shape[0]):
        i_in_var = np.nansum(dmi[cur_sample, :(cur_sample+1)]) + np.nansum(dmi[:(cur_sample+1), cur_sample])
        if i_in_var == 0:
            i_in_var = 1
        i_between_var = np.nansum(dmij[cur_sample, :(cur_sample+1)])
        psrf_like = np.sqrt((i_between_var / (cur_sample+1)) / (i_in_var / (cur_sample+1)))
        psrf_like = [psrf_like]
        for x in range(int(cur_sample * (1 - sample_from)), cur_sample - 1):
            x_in_var = np.nansum(dmi[x, :(cur_sample+1)]) + np.nansum(dmi[:(cur_sample+1), x])
            if x_in_var == 0:
                x_in_var = 1
            x_between_var = np.nansum(dmij[x, :(cur_sample+1)])
            psrf_like.append(np.sqrt((x_between_var / (cur_sample+1)) / (x_in_var / (cur_sample+1))))
        psrf_like = np.median(psrf_like)
        df.append([cur_sample, psrf_like, f"chain_{i}"])

        j_in_var = np.nansum(dmj[cur_sample, :(cur_sample+1)]) + np.nansum(dmj[:(cur_sample+1), cur_sample])
        if j_in_var == 0:
            j_in_var = 1
        j_between_var = np.nansum(dmij[:(cur_sample+1), cur_sample])
        psrf_like = np.sqrt((j_between_var / (cur_sample+1)) / (j_in_var / (cur_sample+1)))
        psrf_like = [psrf_like]
        for x in range(int(cur_sample * (1 - sample_from)), cur_sample - 1):
            x_in_var = np.nansum(dmj[x, :(cur_sample+1)]) + np.nansum(dmj[:(cur_sample+1), x])
            if x_in_var == 0:
                x_in_var = 1
            x_between_var = np.nansum(dmij[:(cur_sample+1), x])
            psrf_like.append(np.sqrt((x_between_var / (cur_sample+1)) / (x_in_var / (cur_sample+1))))
        psrf_like = np.median(psrf_like)
        df.append([cur_sample, psrf_like, f"chain_{j}"])

        if 0.99 < df[-1][1] < 1.01 and 0.99 < df[-2][1] < 1.01:
            consecutive += 1
            for k in cutoff_end.keys():
                if cutoff_end[k] == -1 and consecutive >= int(k * cur_sample):
                    if (cmchain[i].get_pseudo_ess(lower_i=int((cur_sample - consecutive)), upper_i=cur_sample, sample_range=pess_range) >= ess) \
                            and \
                            (cmchain[j].get_pseudo_ess(lower_i=int((cur_sample - consecutive)), upper_i=cur_sample, sample_range=pess_range) >= ess):
                        cutoff_end[k] = cur_sample
                        cutoff_start[k] = int((cur_sample - consecutive))  # * (1 - sample_from))
        else:
            consecutive = 0

    return pd.DataFrame(df, columns=["Sample", "PSRF", "Chain"]), cutoff_start, cutoff_end


def gelman_rubin_trace_plot(cmchain, i, j):
    # Will compute the gelman rubin like diagnostic for chains i and j, this will result in a 'trace'

    sample_from = [0, 0.1, 0.5, 0.75, 1, 0.25, 0.9]
    sample_from = np.sort(sample_from)
    threshold_percentage = [0, 0.1, 0.25, 0.5, 0.75, 1.0]
    threshold_percentage = np.sort(threshold_percentage)

    figure, axis = plt.subplots(ncols=len(sample_from), nrows=len(threshold_percentage),
                                constrained_layout=True,
                                # figsize=[9, 7],
                                squeeze=False, sharex=True, sharey=True)

    for col in range(len(sample_from)):
        df, cutoff_start, cutoff_end = gelman_rubin_trace_with_cutoff(cmchain, i, j, sample_from=sample_from[col],
                                                                      threshold_percentage=threshold_percentage)
        for row in range(len(threshold_percentage)):
            axis[row, col].set_ylim([0.9, 1.1])
            sns.lineplot(data=df, x="Sample", y="PSRF", hue="Chain", alpha=0.5, ax=axis[row, col], legend=False)
            # axis[row, col].axhline(y=1.01, linestyle="--", color="red")
            # axis[row, col].axhline(y=0.99, linestyle="--", color="red")

            axis[row, col].set_xticks = set(df["Sample"])
            # todo the ytick could maybe be adapted
            #  for cases where the value is above, at 1.2 or 1.5 for diverging datasets?
            # axis[row, col].set_ylim([0.8, 1.2])
            # axis[row, col].set_yticks = np.arange(0.8, 1.2, 0.05)

            if cutoff_end[threshold_percentage[row]] != -1 and cutoff_start[threshold_percentage[row]] != -1:
                axis[row, col].axvline(x=cutoff_end[threshold_percentage[row]], color="red")
                # axis[row, col].text(x=cutoff_end[threshold_percentage[row]] + 0.1, y=-.05,
                #                     s=f'{cutoff_end[threshold_percentage[row]]}',
                #                     transform=axis[row, col].get_xaxis_transform(),
                #                     fontsize=8, zorder=20, ha="center",
                #                     va="top", rotation=-45, color="red")
            # if cutoff_start[threshold_percentage[row]] != -1:
                axis[row, col].axvline(x=cutoff_start[threshold_percentage[row]], color="green")
                # axis[row, col].text(x=cutoff_start[threshold_percentage[row]] + 0.1, y=-.5,
                #                     s=f'{cutoff_start[threshold_percentage[row]]}',
                #                     transform=axis[row, col].get_xaxis_transform(),
                #                     fontsize=8, zorder=20, ha="center",
                #                     va="top", rotation=-45, color="green")

                axis[row, col].text(x=cutoff_end[threshold_percentage[row]]-(0.5*(cutoff_end[threshold_percentage[row]] - cutoff_start[threshold_percentage[row]])) + 0.1, y=-.05,
                                    s=f'{cutoff_end[threshold_percentage[row]] - cutoff_start[threshold_percentage[row]] + 1}',
                                    transform=axis[row, col].get_xaxis_transform(),
                                    fontsize=8, zorder=20, ha="center",
                                    va="top", color="black")

            axis[row, col].set_ylabel(f"{threshold_percentage[row]}", color="green")
            axis[row, col].set_xlabel(f"{sample_from[col]}", color="blue")
    figure.supylabel("Threshold Time (fraction of iterations)", color="green")
    figure.supxlabel("Sample from last x-fraction of trees", color="blue")

    # plt.axhline(y=1.05, linestyle="--", color="orange")
    # plt.axhline(y=0.95, linestyle="--", color="orange")
    # 
    # plt.axhline(y=1.1, linestyle="--", color="yellow")
    # plt.axhline(y=0.9, linestyle="--", color="yellow")

    plt.savefig(fname=f"{cmchain.working_dir}/plots/{cmchain.name}_{i}-{j}_grd_singlevalue_evaluation.png",
                format="png", bbox_inches="tight", dpi=800)
    plt.clf()
    plt.close("all")
    gc.collect()


def gr_trace_ess(cmchain, i, j, ess=0, pess_range=100, _overwrite=False, ess_method = "tracerer"):
    # the idea is to set the parameters both to 0.5 which seems to be a good choice in the other evaluations
    tp = 0.5
    sf = 0.5

    df, cutoff_start, cutoff_end = gelman_rubin_trace_with_cutoff(cmchain, i, j, sample_from=sf,
                                                                  threshold_percentage=[tp],
                                                                  ess=ess, pess_range=pess_range)

    figure, axis = plt.subplots(ncols=2, nrows=1,
                                constrained_layout=True, squeeze=False)

    sns.lineplot(data=df, x="Sample", y="PSRF", hue="Chain", alpha=0.5, legend=False, ax=axis[0, 0])

    axis[0, 0].set_xticks = set(df["Sample"])
    axis[0, 0].set_ylim([0.9, 1.1])

    if cutoff_end[tp] != -1 and cutoff_start[tp] != -1:
        axis[0, 0].axvline(x=cutoff_end[tp], color="red")
        axis[0, 0].axvline(x=cutoff_start[tp], color="green")

        axis[0, 0].text(x=cutoff_end[tp] - (
                    0.5 * (cutoff_end[tp] - cutoff_start[tp])) + 0.1, y=.05,
                            s=f'{cutoff_end[tp] - cutoff_start[tp] + 1}',
                            fontsize=8, zorder=20, ha="center",
                            transform=axis[0, 0].get_xaxis_transform(),
                            va="top", color="black")
    axis[0, 0].axhline(y=1.01, linestyle="--", color="red")
    axis[0, 0].axhline(y=0.99, linestyle="--", color="red")

    # Write the cutoff boundaries to a file, if it already exists skip this part
    if _overwrite:
        try:
            os.remove(f"{cmchain.working_dir}/data/{cmchain.name}_{i}_{j}_gress_cutoff{'' if ess == 0 else f'_{ess}'}_{ess_method}")
        except FileNotFoundError:
            pass
    try:
        with open(f"{cmchain.working_dir}/data/{cmchain.name}_{i}_{j}_gress_cutoff{'' if ess == 0 else f'_{ess}'}_{ess_method}", "x") as f:
            f.write(f"{cutoff_start[tp]}\n{cutoff_end[tp]}")
    except FileExistsError:
        pass


    df = _ess_df(cmchain=cmchain, chain_indeces=[i, j], ess_method=ess_method, start=cutoff_start[tp], end=cutoff_end[tp])
    sns.stripplot(data=df, x="Key", y="Value", hue="Chain", ax=axis[0, 1])
    for label in axis[0, 1].get_xticklabels():
        label.set_rotation(90)
    axis[0, 1].set_xlabel("")
    axis[0, 1].set_ylabel("ESS")
    
    axis[0, 1].axhline(y=cutoff_end[tp] - cutoff_start[tp], linestyle="--", color="red")
    axis[0, 1].axhline(y=0.75*(cutoff_end[tp] - cutoff_start[tp]), linestyle="dotted", color="orange")

    plt.savefig(fname=f"{cmchain.working_dir}/plots/{cmchain.name}_{i}-{j}_gress{'' if ess == 0 else f'_{ess}'}_{ess_method}.png",
                format="png", bbox_inches="tight", dpi=400)
    plt.clf()
    plt.close("all")
    gc.collect()

