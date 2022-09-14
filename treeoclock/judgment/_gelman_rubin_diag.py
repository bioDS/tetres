import gc
import os
import random
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from treeoclock.judgment.ess import _ess_df


def _psrf_like_value(dm_in, dm_bt, k, s, e):
    """
        dm_in is the pairwise distance for one tree set - upper triangular form!
        dm_bt is the pairwise distance matrix for two different sets of trees - full matrix
        k is the sample to calculate the psrf value for
        s is the start of the sample of trees
        e is the en of the sample of trees
    """

    in_var = np.nansum(dm_in[k, s:(e+1)]**2) + np.nansum(dm_in[s:(e+1), k]**2)
    bt_var = np.nansum(dm_bt[k, s:(e+1)]**2)

    # todo not sure what to do exactly in this case?
    if in_var == 0:
        return np.sqrt(bt_var/(e-s)+1)
    return np.sqrt(bt_var/in_var)


def gelman_rubin_cut(cmchain, i, j, smoothing, threshold_percentage, ess_threshold=0, pseudo_ess_range=100, ess_method="arviz"):
    # this function will return the cut.start and cut.end values calculated for the given full chain
    dm_i = cmchain.pwd_matrix(i)
    dm_j = cmchain.pwd_matrix(j)
    dm_ij = cmchain.pwd_matrix(i, j)
    dm_ji = np.transpose(dm_ij)

    if dm_i.shape != dm_j.shape:
        raise ValueError("Treesets have different sizes!")

    cutoff_start = -1
    cutoff_end = -1
    consecutive = 0

    for cur_sample in range(dm_i.shape[0]):
        slide_start = int(cur_sample * (1 - smoothing))  # smoothing is impacting the current sliding window start

        psrf_like_i = []
        psrf_like_j = []

        for x in range(slide_start, cur_sample + 1):
            psrf_like_i.append(_psrf_like_value(dm_i, dm_ij, x, slide_start, cur_sample))
            psrf_like_j.append(_psrf_like_value(dm_j, dm_ji, x, slide_start, cur_sample))

        psrf_like_i = np.median(psrf_like_i)
        psrf_like_j = np.median(psrf_like_j)

        if 0.99 < psrf_like_i < 1.01 and 0.99 < psrf_like_j < 1.01:
            consecutive += 1
            if cur_sample - (cur_sample-consecutive) >= ess_threshold:
                if cutoff_end == -1 and consecutive >= int(threshold_percentage * cur_sample):
                        # todo change to calculate the pseudo ess with the already existing distance matrix
                        if (cmchain[i].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample, sample_range=pseudo_ess_range, ess_method=ess_method) >= ess_threshold) \
                                and \
                                (cmchain[j].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample, sample_range=pseudo_ess_range, ess_method=ess_method) >= ess_threshold):
                            cutoff_end = cur_sample
                            cutoff_start = cur_sample - consecutive
                            return cutoff_start, cutoff_end
        else:
            consecutive = 0
    # No cutoff found
    return -1, -1


def gelman_rubin_threshold_list(cmchain, i, j, smoothing, threshold_percentage, ess_threshold=0, pseudo_ess_range=100):
    # This function is able to take a list of threshold_percentage values and also calculates a dataframe of the psrf_like values

    dm_i = cmchain.pwd_matrix(i)
    dm_j = cmchain.pwd_matrix(j)
    dm_ij = cmchain.pwd_matrix(i, j)
    dm_ji = np.transpose(dm_ij)
    
    if dm_i.shape != dm_j.shape:
        raise ValueError("Treesets have different sizes!")

    df = []

    cutoff_start = {k: -1 for k in threshold_percentage}
    cutoff_end = {k: -1 for k in threshold_percentage}
    consecutive = 0

    for cur_sample in range(dm_i.shape[0]):
        slide_start = int(cur_sample * (1 - smoothing))  # smoothing is impacting the current sliding window start
        
        psrf_like_i = []
        psrf_like_j = []
        
        for x in range(slide_start, cur_sample+1):
            psrf_like_i.append(_psrf_like_value(dm_i, dm_ij, x, slide_start, cur_sample))
            psrf_like_j.append(_psrf_like_value(dm_j, dm_ji, x, slide_start, cur_sample))
            
        df.append([cur_sample, np.median(psrf_like_i), f"chain_{i}"])
        df.append([cur_sample, np.median(psrf_like_j), f"chain_{j}"])

        if 0.99 < df[-1][1] < 1.01 and 0.99 < df[-2][1] < 1.01:
            consecutive += 1
            if cur_sample - (cur_sample-consecutive) >= ess_threshold:
                for threshold in cutoff_end.keys():
                    if cutoff_end[threshold] == -1 and consecutive >= int(threshold * cur_sample):
                        # todo change to calculate the pseudo ess with the already existing distance matrix
                        if (cmchain[i].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold) \
                                and \
                                (cmchain[j].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold):
                            cutoff_end[threshold] = cur_sample
                            cutoff_start[threshold] = cur_sample - consecutive
        else:
            consecutive = 0

    return pd.DataFrame(df, columns=["Sample", "PSRF", "Chain"]), cutoff_start, cutoff_end


def gelman_rubin_parameter_choice_plot(cmchain, i, j):
    # Will compute the gelman rubin like diagnostic for chains i and j, this will result in a 'trace'

    # TODO WIP for clean up and rename!
    ess_threshold = 200

    sample_from = [0.1, 0.5, 0.6, 0.75, 1, 0.25, 0.9]
    sample_from = np.sort(sample_from)
    threshold_percentage = [0, 0.1, 0.25, 0.5, 0.75, 1.0]
    threshold_percentage = np.sort(threshold_percentage)

    figure, axis = plt.subplots(ncols=len(sample_from), nrows=len(threshold_percentage),
                                constrained_layout=True,
                                # figsize=[9, 7],
                                squeeze=False, sharex=True, sharey=True)

    for col in range(len(sample_from)):
        df, cutoff_start, cutoff_end = gelman_rubin_threshold_list(cmchain, i, j, smoothing=sample_from[col],
                                                                                threshold_percentage=threshold_percentage,
                                                                                ess_threshold=ess_threshold)
        for row in range(len(threshold_percentage)):
            axis[row, col].set_ylim([0.9, 1.1])
            sns.lineplot(data=df, x="Sample", y="PSRF", hue="Chain", alpha=0.5, ax=axis[row, col], legend=False)

            axis[row, col].set_xticks = set(df["Sample"])

            if cutoff_end[threshold_percentage[row]] != -1 and cutoff_start[threshold_percentage[row]] != -1:
                # adding lines for the cutoff in red and green
                axis[row, col].axvline(x=cutoff_end[threshold_percentage[row]], color="red")
                axis[row, col].axvline(x=cutoff_start[threshold_percentage[row]], color="green")

                # adding text of how many trees are being cut by the specific parameter setting
                axis[row, col].text(x=cutoff_end[threshold_percentage[row]]-(0.5*(cutoff_end[threshold_percentage[row]] - cutoff_start[threshold_percentage[row]])) + 0.1, y=-.05,
                                    s=f'{cutoff_end[threshold_percentage[row]] - cutoff_start[threshold_percentage[row]] + 1}',
                                    transform=axis[row, col].get_xaxis_transform(),
                                    fontsize=8, zorder=20, ha="center",
                                    va="top", color="black")

            axis[row, col].set_ylabel(f"{threshold_percentage[row]}", color="green")
            axis[row, col].set_xlabel(f"{sample_from[col]}", color="blue")

    figure.supylabel("Threshold Time (fraction of iterations)", color="green")
    figure.supxlabel("Sample from last x-fraction of trees", color="blue")

    plt.savefig(fname=f"{cmchain.working_dir}/plots/{cmchain.name}_{i}-{j}_grd_parameter_choices_ess-{ess_threshold}.png",
                format="png", bbox_inches="tight", dpi=800)
    plt.clf()
    plt.close("all")
    gc.collect()


def gelman_rubin_full_chain_subsample(pwts1, pwts2, pwts1ts2, samples: int = 100):
    # this could be a gelman_rubin_like postprocessing tool

    # todo consider burnin as a parameter
    nt1 = pwts1.shape[0]
    nt2 = pwts2.shape[0]
    if not nt1 == nt2:
        raise ValueError("Matrices (Tree sets) should have the same size!")
    psrf_like = []

    if samples != "all":
        for _ in range(samples):
            random_sample = random.randint(0, nt1 - 1)
            r = _psrf_like_value(pwts1, pwts1ts2, random_sample, 0, nt1)
            psrf_like.append(["TS1", r])
            r = _psrf_like_value(pwts2, np.transpose(pwts1ts2), random_sample, 0, nt1)
            psrf_like.append(["TS2", r])
    elif samples == "all":
        for cur_sample in range(nt1):
            r = _psrf_like_value(pwts1, pwts1ts2, cur_sample, 0, nt1)
            psrf_like.append(["TS1", r])
            r = _psrf_like_value(pwts2, np.transpose(pwts1ts2), cur_sample, 0, nt1)
            psrf_like.append(["TS2", r])
    else:
        raise ValueError(f"Unrecognizes samples {samples} given!")
    return pd.DataFrame(data=psrf_like, columns=["Treeset", "PSRF_like"])


def gelman_rubin_all_chains_density_plot(cMChain, samples: int = 100):
    figure, axis = plt.subplots(nrows=cMChain.m_MChains, ncols=cMChain.m_MChains, constrained_layout=True,
                                figsize=[9, 7])

    # initializing the scaled axis for the upper part of the plot:
    scaled_x_axis = np.arange(0.8, 1.7, 0.1)
    scaled_y_axis = np.arange(0, 120, 20)

    for i in range(cMChain.m_MChains - 1):
        for j in range(i + 1, cMChain.m_MChains):
            cur_psrf_like = gelman_rubin_full_chain_subsample(cMChain.pwd_matrix(i),
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
        fname=f"{cMChain.working_dir}/plots/{cMChain.name}_grd_density_full_chain_{'all' if samples == 'all' else f'subsampling_{samples}'}.png",
        format="png", bbox_inches="tight", dpi=800)
    plt.clf()
    plt.close("all")
    gc.collect()
