import gc
import random
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import warnings


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

    if in_var == 0:
        return np.sqrt(bt_var)
    return np.sqrt(bt_var/in_var)


def gelman_rubin_cut(multichain, i, j, smoothing, ess_threshold=200, pseudo_ess_range=100, smoothing_average="mean", _subsampling=False, tolerance=0.02, burnin=0):
    # BEWARE: the values from this function are uncorrected for subsampling, the multichain function is correcting for this
    # this function will return the cut.start and cut.end values calculated for the given full chain

    # The indeces need to be sorted for the dm_ij to make sense, otherwise dm_ij will be the wrong matrix
    if j < i:
        i,j = j, i

    dm_i = multichain.pwd_matrix(i)
    dm_j = multichain.pwd_matrix(j)
    dm_ij = multichain.pwd_matrix(i, j)
    dm_ji = np.transpose(dm_ij)

    # Deleting the given burnin
    dm_i = np.delete(dm_i, list(range(burnin)), axis=0)
    dm_i = np.delete(dm_i, list(range(burnin)), axis=1)
    dm_j = np.delete(dm_j, list(range(burnin)), axis=0)
    dm_j = np.delete(dm_j, list(range(burnin)), axis=1)
    dm_ij = np.delete(dm_ij, list(range(burnin)), axis=0)
    dm_ij = np.delete(dm_ij, list(range(burnin)), axis=1)
    dm_ji = np.delete(dm_ji, list(range(burnin)), axis=0)
    dm_ji = np.delete(dm_ji, list(range(burnin)), axis=1)

    if dm_i.shape != dm_j.shape:
        raise ValueError("Treesets have different sizes!")

    # cutoff_end = -1
    # consecutive = 0
    cutoff_start = 0  # start with the first sample

    if _subsampling:
        for _ in range(_subsampling):
            dm_i = np.delete(dm_i, list(range(0, dm_i.shape[0], 2)), axis=0)
            dm_i = np.delete(dm_i, list(range(0, dm_i.shape[1], 2)), axis=1)

            dm_j = np.delete(dm_j, list(range(0, dm_j.shape[0], 2)), axis=0)
            dm_j = np.delete(dm_j, list(range(0, dm_j.shape[1], 2)), axis=1)

            dm_ij = np.delete(dm_ij, list(range(0, dm_ij.shape[0], 2)), axis=0)
            dm_ij = np.delete(dm_ij, list(range(0, dm_ij.shape[1], 2)), axis=1)

            dm_ji = np.delete(dm_ji, list(range(0, dm_ji.shape[0], 2)), axis=0)
            dm_ji = np.delete(dm_ji, list(range(0, dm_ji.shape[1], 2)), axis=1)

    for cur_sample in range(dm_i.shape[0]):
        slide_start = int(cur_sample * (1 - smoothing))  # smoothing is impacting the current sliding window start

        # psrf_like_i = _psrf_like_value(dm_i, dm_ij, k=cur_sample, s=cutoff_start, e=cur_sample)
        # psrf_like_j = _psrf_like_value(dm_j, dm_ji, k=cur_sample, s=cutoff_start, e=cur_sample)

        psrf_like_i = []
        psrf_like_j = []
        for x in range(slide_start, cur_sample + 1):
            psrf_like_i.append(_psrf_like_value(dm_i, dm_ij, x, slide_start, cur_sample))
            psrf_like_j.append(_psrf_like_value(dm_j, dm_ji, x, slide_start, cur_sample))
        if smoothing_average == "median":
            psrf_like_i = np.median(psrf_like_i)
            psrf_like_j = np.median(psrf_like_j)
        elif smoothing_average == "mean":
            psrf_like_i = np.mean(psrf_like_i)
            psrf_like_j = np.mean(psrf_like_j)
        else:
            raise ValueError(f"Unrecognized parameter {smoothing_average} smoothing_average!")
        # print(psrf_like_i, psrf_like_j)
        if 1/(1 + tolerance) < psrf_like_i < 1+tolerance and 1/(1 + tolerance) < psrf_like_j < 1+tolerance:
            if cur_sample - cutoff_start > ess_threshold:
                    # todo change to calculate the pseudo ess with the already existing distance matrix
                    if (multichain[i].get_pseudo_ess(lower_i=cutoff_start, upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold) \
                            and \
                            (multichain[j].get_pseudo_ess(lower_i=cutoff_start, upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold):
                    # if (multichain[i].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold) \
                    #         and \
                    #         (multichain[j].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold):
                        # cutoff_end = cur_sample
                        # cutoff_start = cur_sample - consecutive
                        return cutoff_start+burnin, cur_sample+burnin
        else:
            # consecutive = 0
            # if we are outside the boundary we discard the previous samples from the cut as they don't count
            cutoff_start = cur_sample
    # No cutoff found
    return -1, -1


def gelman_rubin_ess_threshold_list_list(multichain, i, j, smoothing, ess_threshold_list, pseudo_ess_range=100, smoothing_average="median", _subsampling=False, tolerance=0.02):
    # This function is able to take a list of threshold_percentage values and also calculates a dataframe of the psrf_like values

    warnings.warn("This funciton has not been updated to the latest version of the GRdiagnostic yet")

    dm_i = multichain.pwd_matrix(i)
    dm_j = multichain.pwd_matrix(j)
    dm_ij = multichain.pwd_matrix(i, j)
    dm_ji = np.transpose(dm_ij)
    
    if dm_i.shape != dm_j.shape:
        raise ValueError("Treesets have different sizes!")

    df = []

    cutoff_start = {k: 0 for k in ess_threshold_list}
    cutoff_end = {k: -1 for k in ess_threshold_list}
    # consecutive = 0

    if _subsampling:
        # delete every second row (can be repeated multiple times)
        for _ in range(_subsampling):
            dm_i = np.delete(dm_i, list(range(0, dm_i.shape[0], 2)), axis=0)
            dm_i = np.delete(dm_i, list(range(0, dm_i.shape[1], 2)), axis=1)

            dm_j = np.delete(dm_j, list(range(0, dm_j.shape[0], 2)), axis=0)
            dm_j = np.delete(dm_j, list(range(0, dm_j.shape[1], 2)), axis=1)

            dm_ij = np.delete(dm_ij, list(range(0, dm_ij.shape[0], 2)), axis=0)
            dm_ij = np.delete(dm_ij, list(range(0, dm_ij.shape[1], 2)), axis=1)

            dm_ji = np.delete(dm_ji, list(range(0, dm_ji.shape[0], 2)), axis=0)
            dm_ji = np.delete(dm_ji, list(range(0, dm_ji.shape[1], 2)), axis=1)

    for cur_sample in range(0, dm_i.shape[0]):
        slide_start = int(cur_sample * (1 - smoothing))  # smoothing is impacting the current sliding window start
        
        psrf_like_i = []
        psrf_like_j = []
        
        for x in range(slide_start, cur_sample+1):
            psrf_like_i.append(_psrf_like_value(dm_i, dm_ij, x, slide_start, cur_sample))
            psrf_like_j.append(_psrf_like_value(dm_j, dm_ji, x, slide_start, cur_sample))

        if smoothing_average == "median":
            df.append([cur_sample, np.median(psrf_like_i), f"chain_{i}"])
            df.append([cur_sample, np.median(psrf_like_j), f"chain_{j}"])
        elif smoothing_average == "mean":
            df.append([cur_sample, np.mean(psrf_like_i), f"chain_{i}"])
            df.append([cur_sample, np.mean(psrf_like_j), f"chain_{j}"])
        else:
            raise ValueError(f"Smoothing_function = {smoothing_average} not recognized!")

        if 1/(1 + tolerance) < df[-1][1] < 1+tolerance and 1/(1 + tolerance) < df[-2][1] < 1+tolerance:
            # consecutive += 1
            for ess_threshold in cutoff_end.keys():
                if cutoff_start[ess_threshold] == -1:
                    cutoff_start[ess_threshold] = slide_start
                # if cutoff_end[ess_threshold] == -1 and consecutive >= ess_threshold:

                if cutoff_end[ess_threshold] == -1 and cur_sample - cutoff_start[ess_threshold] > ess_threshold:
                    # todo change to calculate the pseudo ess with the already existing distance matrix
                    if (multichain[i].get_pseudo_ess(lower_i=cutoff_start[ess_threshold], upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold) \
                            and \
                            (multichain[j].get_pseudo_ess(lower_i=cutoff_start[ess_threshold], upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold):
                        cutoff_end[ess_threshold] = cur_sample
                        # cutoff_start[ess_threshold] = cur_sample - consecutive
        else:
            for ess_threshold in cutoff_end.keys():
                if cutoff_end[ess_threshold] == -1:
                    cutoff_start[ess_threshold] = -1
            # consecutive = 0

    # todo without the consecutive stuff, it can happen that cutoffstart is set but the end is never found
    #  renormalizing to fit previous strucutre where this was not possible
    # for ess_threshold in cutoff_end.keys():
    #     if cutoff_end[ess_threshold] == -1:
    #         cutoff_start[ess_threshold] = -1

    return pd.DataFrame(df, columns=["Sample", "PSRF", "Chain"]), cutoff_start, cutoff_end


def gelman_rubin_parameter_choice_plot(multichain, i, j, _subsampling=False, tolerance=0.02, smoothing_average="mean"):
    warnings.warn("This funciton has not been updated to the latest version of the GRdiagnostic yet")

    warnings.warn("This takes a long time and is not optimized, nor does it save anything other than the plot!")
    # Will compute the gelman rubin like diagnostic for chains i and j, this will result in a 'trace'
    ess_threshold_list = np.sort([200, 500])
    smoothing = [0.6, 0.9, 1]
    smoothing = np.sort(smoothing)

    figure, axis = plt.subplots(ncols=len(smoothing), nrows=len(ess_threshold_list),
                                constrained_layout=True,
                                # figsize=[9, 7],
                                squeeze=False, sharex=True, sharey=True)

    for col in range(len(smoothing)):
        df, cutoff_start, cutoff_end = gelman_rubin_ess_threshold_list_list(multichain, i, j, smoothing=smoothing[col],
                                                                            ess_threshold_list=ess_threshold_list,
                                                                            smoothing_average=smoothing_average,
                                                                            _subsampling=_subsampling,
                                                                            tolerance=tolerance)
        for row in range(len(ess_threshold_list)):
            axis[row, col].set_ylim([0.9, 1.1])
            sns.lineplot(data=df, x="Sample", y="PSRF", hue="Chain", alpha=0.5, ax=axis[row, col], legend=False, linewidth=0.6)

            axis[row, col].set_xticks = set(df["Sample"])

            if cutoff_end[ess_threshold_list[row]] != -1 and cutoff_start[ess_threshold_list[row]] != -1:
                # adding lines for the cutoff in red and green
                axis[row, col].axvline(x=cutoff_end[ess_threshold_list[row]], color="red", linewidth=0.8)
                axis[row, col].axvline(x=cutoff_start[ess_threshold_list[row]], color="green", linewidth=0.8)

                # adding text of how many trees are being cut by the specific parameter setting
                axis[row, col].text(x=cutoff_end[ess_threshold_list[row]]-(0.5*(cutoff_end[ess_threshold_list[row]] - cutoff_start[ess_threshold_list[row]])) + 0.1, y=-.05,
                                    s=f'{cutoff_end[ess_threshold_list[row]] - cutoff_start[ess_threshold_list[row]] + 1}',
                                    transform=axis[row, col].get_xaxis_transform(),
                                    fontsize=8, zorder=20, ha="center",
                                    va="top", color="black")

            axis[row, col].set_ylabel(f"{ess_threshold_list[row]}", color="green")
            axis[row, col].set_xlabel(f"{smoothing[col]}", color="blue")

    figure.supylabel("Pseudo ESS threshold", color="green")
    figure.supxlabel("Smoothing fraction", color="blue")

    plt.savefig(fname=f"{multichain.working_dir}/plots/{multichain.name}_{i}-{j}_grd_parameter_choices{f'_subsampling-{_subsampling}' if _subsampling else ''}_ess-list_{smoothing_average}_{tolerance}.pdf",
                format="pdf", bbox_inches="tight", dpi=1200)
    plt.clf()
    plt.close("all")
    gc.collect()


def gelman_rubin_full_chain_subsample(pwts1, pwts2, pwts1ts2, samples: int = 100, no_smooth=False):
    # this could be a gelman_rubin_like postprocessing tool
    nt1 = pwts1.shape[0]
    nt2 = pwts2.shape[0]
    if not nt1 == nt2:
        raise ValueError("Matrices (Tree sets) should have the same size!")
    psrf_like = []

    if samples != "all":
        for x in range(samples):  # todo sample outside to have proper sample indeces for the dataframe....
            random_sample = random.randint(0, nt1 - 1)
            r = _psrf_like_value(pwts1, pwts1ts2, random_sample, 0, nt1)
            psrf_like.append(["TS1", r, x])
            r = _psrf_like_value(pwts2, np.transpose(pwts1ts2), random_sample, 0, nt1)
            psrf_like.append(["TS2", r, x])
    elif samples == "all":
        s = 0
        e = nt1
        for cur_sample in range(nt1):
            if no_smooth:
                s = cur_sample
                e = cur_sample
            r = _psrf_like_value(pwts1, pwts1ts2, cur_sample, s, e)
            psrf_like.append(["TS1", r, cur_sample])
            r = _psrf_like_value(pwts2, np.transpose(pwts1ts2), cur_sample, s, e)
            psrf_like.append(["TS2", r, cur_sample])
    else:
        raise ValueError(f"Unrecognizes samples {samples} given!")
    return pd.DataFrame(data=psrf_like, columns=["Treeset", "PSRF_like", "Sample"])


def gelman_rubin_all_chains_density_plot(multichain, samples: int = 100):
    figure, axis = plt.subplots(nrows=multichain.m_chains, ncols=multichain.m_chains, layout="constrained",
                                figsize=[9, 7])

    # Diagonal Plots deleted as they are not used
    for i in range(multichain.m_chains):
        figure.delaxes(axis[i, i])

    # initializing the scaled axis for the upper part of the plot:
    scaled_x_axis = np.arange(0.8, 1.7, 0.1)
    scaled_y_axis = np.arange(0, 120, 20)

    for i in range(multichain.m_chains - 1):
        for j in range(i + 1, multichain.m_chains):
            cur_psrf_like = gelman_rubin_full_chain_subsample(multichain.pwd_matrix(i),
                                                              multichain.pwd_matrix(j),
                                                              multichain.pwd_matrix(i, j),
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

    # Annotation of the plot
    pad = 5  # in points
    cols = [f"chain{m}" for m in range(multichain.m_chains)]
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
        fname=f"{multichain.working_dir}/plots/{multichain.name}_grd_density_full_chain_{'all' if samples == 'all' else f'subsampling_{samples}'}.pdf",
        format="pdf", bbox_inches="tight", dpi=1200)
    plt.clf()
    plt.close("all")
    gc.collect()


def density_trace_plot(multichain, interval, i=0, j=1, no_smooth=False):

    # figure, axis = plt.subplots(nrows=1, ncols=2, layout="constrained", figsize=(15, 5))
    figure, axis = plt.subplot_mosaic([["a", "a", "b"]], figsize=(15, 5), sharey=True)
    ax = [(label, ax) for label, ax in axis.items()]

    cur_psrf_like = gelman_rubin_full_chain_subsample(multichain.pwd_matrix(i),
                                                              multichain.pwd_matrix(j),
                                                              multichain.pwd_matrix(i, j),
                                                              samples="all", no_smooth=no_smooth)

    sns.kdeplot(ax=ax[1][1],
                data=cur_psrf_like,
                y="PSRF_like",
                hue="Treeset",
                fill=True,
                legend=False,
                common_norm=False)

    x_ticks = [(i+1)*interval for i in list(cur_psrf_like["Sample"])]
    sns.lineplot(ax=ax[0][1], data=cur_psrf_like, y="PSRF_like", x=x_ticks, hue="Treeset")

    ax[0][1].set_xlim((0, max(x_ticks)))

    new_ticks = np.linspace(0, max(x_ticks), 11, endpoint=True, dtype=int)
    ax[0][1].set_xticks(new_ticks, new_ticks, fontsize=12, rotation=45)

    ax[0][1].set_xlabel("Sample", fontsize=12)
    ax[1][1].set_xlabel("Density", fontsize=12)
    ax[0][1].grid()
    ax[1][1].grid()

    ax[1][1].tick_params(axis="y", which="both", left=False)

    ax[0][1].set_ylabel("PSRF", fontsize=12)
    ax[0][1].legend(ncol=2)

    ax[0][1].tick_params(axis="both", labelsize=12)
    ax[1][1].tick_params(axis="both", labelsize=12)

    # plt.subplots_adjust(wspace=0.01)
    plt.tight_layout()
    # plt.show()
    plt.savefig(
        fname=f"{multichain.working_dir}/plots/{multichain.name}_density_trace_plot{'_no-smooth' if no_smooth else ''}_{i}-{j}.pdf",
        format="pdf", bbox_inches="tight", dpi=1200)
    plt.clf()
    plt.close("all")
    gc.collect()
