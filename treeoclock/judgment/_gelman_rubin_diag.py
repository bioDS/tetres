from treeoclock.trees.time_trees import TimeTreeSet, findpath_distance
from treeoclock.judgment._pairwise_distance_matrix import calc_pw_distances, calc_pw_distances_two_sets

import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def gelman_rubin_distance_diagnostic_plot(cMChain, samples: int = 100):

    figure, axis = plt.subplots(nrows=cMChain.m_MChains, ncols=cMChain.m_MChains, constrained_layout=True, figsize=[9, 7])
    for i in range(cMChain.m_MChains-1):
        for j in range(i+1, cMChain.m_MChains):
            cur_psrf_like = gelman_rubin_distance_diagnostic_from_matrices(cMChain.pwd_matrix(i),
                                                                           cMChain.pwd_matrix(j),
                                                                           cMChain.pwd_matrix(i, j),
                                                                           samples=samples)
            sns.kdeplot(ax=axis[i, j], data=cur_psrf_like, fill=True)
            for label in axis[i, j].get_xticklabels():
                label.set_rotation(45)
            sns.boxplot(ax=axis[j, i], data=cur_psrf_like)
            # todo axis[i, i] some plot on diagonal missing

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

    plt.savefig(fname=f"{cMChain.working_dir}/{cMChain.name}_grd_plot.png", format="png", bbox_inches="tight", dpi=800)


def gelman_rubin_distance_diagnostic_from_matrices(pwts1, pwts2, pwts1ts2,
                                     samples: int = 100):
    nt1 = pwts1.shape[0]
    nt2 = pwts2.shape[0]
    psrf_like = []

    for _ in range(samples):
        ts1_sample = random.randint(0, nt1 - 1)
        in_sample_var = np.sum(pwts1[ts1_sample, :]) + np.sum(pwts1[:, ts1_sample])
        between_sample_var = np.sum(pwts1ts2[ts1_sample, :])
        psrf_like.append(np.sqrt((between_sample_var / nt2) / (in_sample_var / nt1)))

        ts2_sample = random.randint(0, nt2 - 1)
        in_sample_var = np.sum(pwts2[ts2_sample, :]) + np.sum(pwts2[:, ts2_sample])
        between_sample_var = np.sum(pwts1ts2[:, ts2_sample])
        psrf_like.append(np.sqrt((between_sample_var / nt1) / (in_sample_var / nt2)))

    return psrf_like


def gelman_rubin_distance_diagnostic(tree_set1: TimeTreeSet, tree_set2: TimeTreeSet,
                                     samples: int = 100):
    # todo time the distance matrix computation and text remco how long it takes to do this
    distances_pwts1 = calc_pw_distances(tree_set1)
    distances_pwts2 = calc_pw_distances(tree_set2)
    distance_pwts1ts2 = calc_pw_distances_two_sets(tree_set1, tree_set2)

    # var_difference = []
    # in_sample_var_list = []
    # between_sample_var_list = []
    psrf_like = []

    for _ in range(samples):
        ts1_sample = random.randint(0, len(tree_set1) - 1)

        in_sample_var = np.sum(distances_pwts1[ts1_sample, :]) + np.sum(distances_pwts1[:, ts1_sample])

        between_sample_var = np.sum(distance_pwts1ts2[ts1_sample, :])

        # in_sample_var_list.append(in_sample_var / len(tree_set1))
        # between_sample_var_list.append(between_sample_var / len(tree_set2))
        # var_difference.append(abs((between_sample_var / len(tree_set2)) - (in_sample_var / len(tree_set1))))
        psrf_like.append(np.sqrt((between_sample_var / len(tree_set2)) / (in_sample_var / len(tree_set1))))

        # Same for a sample from the second tree set
        ts2_sample = random.randint(0, len(tree_set2) - 1)
        in_sample_var = np.sum(distances_pwts2[ts2_sample, :]) + np.sum(distances_pwts2[:, ts2_sample])
        between_sample_var = np.sum(distance_pwts1ts2[:, ts2_sample])
        psrf_like.append(np.sqrt((between_sample_var / len(tree_set1)) / (in_sample_var / len(tree_set2))))

    # return var_difference, in_sample_var_list, between_sample_var_list, psrf_like
    return psrf_like