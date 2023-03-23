import gc

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

global _gr_boundary
_gr_boundary = 0.02


def _psrf_like_fast_value(dm_in, dm_bt, k, s, e, fix_trees, weights_in, weights_other):
    """
        dm_in is the pairwise distance for one tree set - upper triangular form!
        dm_bt is the pairwise distance matrix for two different sets of trees - full matrix
        k is the sample to calculate the psrf value for
        s is the start of the sample of trees
        e is the end of the sample of trees
    """

    in_var = 0
    bt_var = 0

    for f in fix_trees:
        if s <= f <= e:
            # in_var = np.nansum(dm_in[k, s:(e+1)]**2) + np.nansum(dm_in[s:(e+1), k]**2)
            in_var += ((dm_in[k, f] + dm_in[f, k]) * weights_in[f])**2

            # bt_var = np.nansum(dm_bt[k, s:(e+1)]**2)
            bt_var += (dm_bt[k, f] * weights_other[f])**2

    if in_var == 0:
        return np.sqrt(bt_var)
    return np.sqrt(bt_var/in_var)


def fast_grd_cut(cmchain, i, j, smoothing, ess_threshold=0, pseudo_ess_range=100, smoothing_average="median",
                 fix_tree_fraction=0.8):

    if len(cmchain[i]) != len(cmchain[j]):
        raise ValueError("Tree sets need to have the same size!")

    # this function will return the cut.start and cut.end values calculated for the given full chain
    global _gr_boundary


    # Creating 100 fixed points over the length of treesets i and j
    #  in the final implementation this 100 should instead be computed on the go and not be preset!
    if len(cmchain[i]) > 1001:
        print(
            "Current Implementation is for 1000 trees in the set, "
            "if there are more there needs to be adjusting to the fixed tree indeces!")

    # fix_point_fraction = 1  # Number of fixed trees as a fraction of the full length
    fix_points = int(len(cmchain[i].trees) * fix_tree_fraction)

    fix_trees = np.linspace(start=0, stop=len(cmchain[i]) - 1, num=int(fix_points), dtype=int)
    weights_i = {t: 1.0 for t in fix_trees}
    weights_j = {t: 1.0 for t in fix_trees}

    dm_i = cmchain.pwd_matrix(i)
    dm_j = cmchain.pwd_matrix(j)
    dm_ij = cmchain.pwd_matrix(i, j)
    dm_ji = np.transpose(dm_ij)

    m = 4  # fuzzy cluster m, larger values make the clusteres fuzzier

    ###
    # Inefficient precomputation of the weights, but thats not the point here!
    # Todo should be its own funciton at some point, also use for centroid
    ###

    # Compute the weights for trees
    for cur_sample in range(dm_i.shape[0]):
        slide_start = int(cur_sample * (1 - smoothing))  # smoothing is impacting the current sliding window start
        # slide window is from slide_start to slide end!

        if cur_sample in fix_trees:
            # weights are initialized with 1, so this can be skipped
            continue

        upto_now_fixed = fix_trees[(slide_start <= fix_trees) & (fix_trees < cur_sample)]
        cur_dist_i = {}
        cur_dist_j = {}

        i_zero = False
        j_zero = False

        for update_tree in upto_now_fixed:
            if dm_i[update_tree, cur_sample] == 0:
                i_zero = True
                weights_i[update_tree] += 1
            elif not i_zero:
                cur_dist_i[update_tree] = dm_i[update_tree, cur_sample]  # + dm_j[cur_sample, update_fix]

            if dm_j[update_tree, cur_sample] == 0:
                j_zero = True
                weights_j[update_tree] += 1
            elif not j_zero:
                cur_dist_j[update_tree] = dm_j[update_tree, cur_sample]  # + dm_j[cur_sample, update_fix]

        if i_zero and j_zero:
            # this is for cases where in both chains the same tree as the respective fixed tree got sampled
            continue

        # s_i = sum(cur_dist_i.values())
        # s_j = sum(cur_dist_j.values())

        for update_tree in upto_now_fixed:
            # tmp_sum_i = (1 - prev_prop_i)
            # tmp_sum_j = (1 - prev_prop_j) - (cur_dist_j[update_fix] / s_j)
            if not i_zero:
                weights_i[update_tree] += (1 / sum(
                    [((cur_dist_i[update_tree] / cur_dist_i[f]) ** (2 / (m - 1))) for f in
                     upto_now_fixed]))  # (cur_dist_i[update_fix] / s_i) ** (2/m-1))
            if not j_zero:
                weights_j[update_tree] += (1 / sum(
                    [((cur_dist_j[update_tree] / cur_dist_j[f]) ** (2 / (m - 1))) for f in upto_now_fixed]))

    consecutive = 0

    for cur_sample in range(dm_i.shape[0]):
        slide_start = int(cur_sample * (1 - smoothing))  # smoothing is impacting the current sliding window start

        psrf_like_i = []
        psrf_like_j = []

        for x in range(slide_start, cur_sample + 1):
            psrf_like_i.append(_psrf_like_fast_value(dm_i, dm_ij, x, slide_start, cur_sample, fix_trees, weights_i, weights_j))
            psrf_like_j.append(_psrf_like_fast_value(dm_j, dm_ji, x, slide_start, cur_sample, fix_trees, weights_j, weights_i))

        if smoothing_average == "median":
            psrf_like_i = np.median(psrf_like_i)
            psrf_like_j = np.median(psrf_like_j)
        elif smoothing_average == "mean":
            psrf_like_i = np.mean(psrf_like_i)
            psrf_like_j = np.mean(psrf_like_j)
        else:
            raise ValueError(f"Unrecognized parameter {smoothing_average} smoothing_average!")

        if 1-_gr_boundary < psrf_like_i < 1+_gr_boundary and 1-_gr_boundary < psrf_like_j < 1+_gr_boundary:
            consecutive += 1
            if cur_sample - (cur_sample - consecutive) >= ess_threshold and (
                    cmchain[i].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample,
                                              sample_range=pseudo_ess_range) >= ess_threshold) \
                    and \
                    (cmchain[j].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample,
                                               sample_range=pseudo_ess_range) >= ess_threshold):
                # if cutoff_end == -1 and consecutive >= int(threshold_percentage * cur_sample):
                        # todo change to calculate the pseudo ess with the already existing distance matrix
                        # if (cmchain[i].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold) \
                        #         and \
                        #         (cmchain[j].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample, sample_range=pseudo_ess_range) >= ess_threshold):
                    cutoff_end = cur_sample
                    cutoff_start = cur_sample - consecutive
                    return cutoff_start, cutoff_end
        else:
            consecutive = 0
    # No cutoff found
    return -1, -1



def fast_approx_fix_fraction_plot(cmchain, i, j):

    # todo add the original diagnostic cut also to double check that ftf=1 is the same as this!

    smt = 0.5
    esst = 200
    pessr = 100
    smtavg = "mean"

    # todo multiprocess this part i guess?
    for ftf in np.arange(0.1, 1.1, 0.1):

        cur_start, cur_end = fast_grd_cut(cmchain, i, j, smoothing=smt, ess_threshold=esst,
                              pseudo_ess_range=pessr, smoothing_average=smtavg, fix_tree_fraction=ftf)
        print(ftf, cur_start, cur_end)


def gelman_rubin_fast_approx_threshold_list(cmchain, i, j, smoothing, ess_threshold=0, pseudo_ess_range=100,
                                smoothing_average="median"):

    if len(cmchain[i]) != len(cmchain[j]):
        raise ValueError("Tree sets need to have the same size!")

    # This function is able to take a list of threshold_percentage values and also calculates a dataframe of the psrf_like values
    global _gr_boundary

    # precompute distance matrices
    dm_i = cmchain.pwd_matrix(i)
    dm_j = cmchain.pwd_matrix(j)
    dm_ij = cmchain.pwd_matrix(i, j)
    dm_ji = np.transpose(dm_ij)

    # todo investigate this a little more
    m = 4  # fuzzy cluster m, larger values make the clusteres fuzzier

    fix_point_fraction = 1  # Number of fixed trees as a fraction of the full length
    fix_points = int(len(cmchain[i].trees) * fix_point_fraction)

    # todo if any of the fixed trees have distance 0 then one should be eliminated and switched for something else!
    # fix_trees = np.linspace(start=0, stop=len(cmchain[i])-1,num = fix_points, dtype=int)
    fix_trees = np.concatenate((np.linspace(start=0, stop=int((len(cmchain[i]) - 1) / 2) - 1, num=int(fix_points * 0.3), dtype=int),
                                np.linspace(start=int((len(cmchain[i]) - 1) / 2), stop=len(cmchain[i]) - 1, num=int(fix_points * 0.7), dtype=int)))

    # delete duplicate fix trees, i.e. dm[i,j]=0
    change_i = set()
    change_j = set()
    for index_i, ei in enumerate(fix_trees):
        for index_j, ej in enumerate(fix_trees[i+1:]):
            if dm_i[ei, ej] == 0:
                change_i.add(index_j+i)
            if dm_j[ei, ej] == 0:
                change_j.add(index_j+i)

    if len(change_i) != 0:
        print('Duplicated fixes!')  # todo
    if len(change_j) != 0:
        print('Duplicated fixes!')  # todo

    # todo investiage how the number impacts the result of the heuristic
    # todo there could be some threshold at which a new tree is added as a fixed tree instead of weighting it!

    weights_i = {t: 1.0 for t in fix_trees}
    weights_j = {t: 1.0 for t in fix_trees}


    ###
    # Inefficient precomputation of the weights, but thats not the point here!
    # Todo should be its own funciton at some point, also use for centroid
    ###

    # Compute the weights for trees
    for cur_sample in range(dm_i.shape[0]):
        slide_start = int(cur_sample * (1 - smoothing))  # smoothing is impacting the current sliding window start
        # slide window is from slide_start to slide end!

        if cur_sample in fix_trees:
            # weights are initialized with 1, so this can be skipped
            continue

        upto_now_fixed = fix_trees[(slide_start <= fix_trees) & (fix_trees < cur_sample)]
        cur_dist_i = {}
        cur_dist_j = {}

        i_zero = False
        j_zero = False

        # todo speed this function up

        for update_tree in upto_now_fixed:
            if dm_i[update_tree, cur_sample] == 0:
                i_zero = True
                weights_i[update_tree] += 1
            elif not i_zero:
                cur_dist_i[update_tree] = dm_i[update_tree, cur_sample]  # + dm_j[cur_sample, update_fix]

            if dm_j[update_tree, cur_sample] == 0:
                j_zero = True
                weights_j[update_tree] += 1
            elif not j_zero:
                cur_dist_j[update_tree] = dm_j[update_tree, cur_sample]  # + dm_j[cur_sample, update_fix]

        if i_zero and j_zero:
            # this is for cases where in both chains the same tree as the respective fixed tree got sampled
            continue

        # s_i = sum(cur_dist_i.values())
        # s_j = sum(cur_dist_j.values())

        for update_tree in upto_now_fixed:
            # tmp_sum_i = (1 - prev_prop_i)
            # tmp_sum_j = (1 - prev_prop_j) - (cur_dist_j[update_fix] / s_j)
            if not i_zero:
                weights_i[update_tree] += (1 / sum([((cur_dist_i[update_tree] /cur_dist_i[f]) ** (2/(m-1))) for f in upto_now_fixed]))  # (cur_dist_i[update_fix] / s_i) ** (2/m-1))
            if not j_zero:
                weights_j[update_tree] += (1 / sum([((cur_dist_j[update_tree] /cur_dist_j[f]) ** (2/(m-1))) for f in upto_now_fixed]))

    df = []

    cutoff_start = -1
    cutoff_end = -1
    consecutive = 0

    for cur_sample in range(dm_i.shape[0]):
        slide_start = int(cur_sample * (1 - smoothing))  # smoothing is impacting the current sliding window start

        psrf_like_i = []
        psrf_like_j = []

        for x in range(slide_start, cur_sample + 1):
            psrf_like_i.append(_psrf_like_fast_value(dm_i, dm_ij, x, slide_start, cur_sample, fix_trees, weights_i, weights_j))
            psrf_like_j.append(_psrf_like_fast_value(dm_j, dm_ji, x, slide_start, cur_sample, fix_trees, weights_j, weights_i))

        if smoothing_average == "median":
            df.append([cur_sample, np.median(psrf_like_i), f"chain_{i}"])
            df.append([cur_sample, np.median(psrf_like_j), f"chain_{j}"])
        elif smoothing_average == "mean":
            df.append([cur_sample, np.mean(psrf_like_i), f"chain_{i}"])
            df.append([cur_sample, np.mean(psrf_like_j), f"chain_{j}"])
        else:
            raise ValueError(f"Smoothing_function = {smoothing_average} not recognized!")

        if 1/(1+_gr_boundary) < df[-1][1] < 1 + _gr_boundary and 1/(1+_gr_boundary) < df[-2][1] < 1 + _gr_boundary:
            consecutive += 1
            if cur_sample - (cur_sample - consecutive) >= ess_threshold:
                if cutoff_end == -1:
                    # todo change to calculate the pseudo ess with the already existing distance matrix
                    if (cmchain[i].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample,
                                                    sample_range=pseudo_ess_range) >= ess_threshold) \
                            and \
                            (cmchain[j].get_pseudo_ess(lower_i=cur_sample - consecutive, upper_i=cur_sample,
                                                        sample_range=pseudo_ess_range) >= ess_threshold):
                        cutoff_end = cur_sample
                        cutoff_start = cur_sample - consecutive
        else:
            consecutive = 0

    return pd.DataFrame(df, columns=["Sample", "PSRF", "Chain"]), cutoff_start, cutoff_end


def gelman_rubin_fast_approx_parameter_choice_plot(cmchain, i, j):
    # Will compute the gelman rubin like diagnostic for chains i and j, this will result in a 'trace'
    # TODO WIP for clean up and rename!
    ess_threshold = 200
    smoothing_function = "mean"  # TODO

    # sample_from = [0.1, 0.5, 0.6, 0.75, 1, 0.25, 0.9]
    sample_from = [0.5, 0.75]
    sample_from = np.sort(sample_from)

    figure, axis = plt.subplots(ncols=len(sample_from), nrows=1,
                                constrained_layout=True, sharey=True)
                                # figsize=[9, 7],
                                # squeeze=False)

    for col in range(len(sample_from)):
        df, cutoff_start, cutoff_end = gelman_rubin_fast_approx_threshold_list(cmchain, i, j, smoothing=sample_from[col],
                                                                   ess_threshold=ess_threshold,
                                                                   smoothing_average=smoothing_function)
        # todo set the ylim for the plot
        # axis.set_ylim([0.9, 1.1])  # todo this doesnt work....
        sns.lineplot(data=df, x="Sample", y="PSRF", hue="Chain", alpha=0.5, ax=axis[col], legend=False)

        # axis[col].set_xticks = set(df["Sample"])

        if cutoff_end != -1 and cutoff_start != -1:
            # adding lines for the cutoff in red and green
            axis[col].axvline(x=cutoff_end, color="red")
            axis[col].axvline(x=cutoff_start, color="green")

            # adding text of how many trees are being cut by the specific parameter setting
            # axis[col].text(x=cutoff_end - (0.5 * (
            #             cutoff_end - cutoff_start)) + 0.1,
            #                     y=-.05,
            #                     s=f'{cutoff_end - cutoff_start + 1}',
            #                     transform=axis[col].get_xaxis_transform(),
            #                     fontsize=8, zorder=20, ha="center",
            #                     va="top", color="black")

        # axis[col].set_ylabel(f"{[row]}", color="green")
        axis[col].set_xlabel(f"{sample_from[col]}", color="blue")

    figure.supylabel("Threshold Time (fraction of iterations)", color="green")
    figure.supxlabel("Sample from last x-fraction of trees", color="blue")

    plt.savefig(
        fname=f"{cmchain.working_dir}/plots/{cmchain.name}_{i}-{j}_grd_fast_approximation_parameter_choices_ess-{ess_threshold}_{smoothing_function}.png",
        format="png", bbox_inches="tight", dpi=800)
    plt.clf()
    plt.close("all")
    gc.collect()
