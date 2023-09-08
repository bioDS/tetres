import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tetres.trees.time_trees import TimeTreeSet


def discrete_cladeset_comparator(tree_set_i, tree_set_j, plot=True, burnin=0, file=""):
    internal_i, internal_j = TimeTreeSet(), TimeTreeSet()

    internal_i.map = tree_set_i.map
    internal_j.map = tree_set_j.map

    internal_i.trees = tree_set_i.trees[int(len(tree_set_i)*burnin):]
    internal_j.trees = tree_set_j.trees[int(len(tree_set_j)*burnin):]

    i_clade_dict = internal_i.get_clade_rank_dictionary()
    j_clade_dict = internal_j.get_clade_rank_dictionary()

    joint_keys = i_clade_dict.keys() & j_clade_dict.keys()  # overlap of keys
    all_keys = i_clade_dict.keys() | j_clade_dict.keys()  # combined set of keys

    df = []
    line_coords = []

    count_outsides = 0
    for k in all_keys:
        # compute the percentage within the diagonals

        if k in joint_keys:
            x = np.mean(i_clade_dict[k])
            y = np.mean(j_clade_dict[k])
            s = np.mean([len(i_clade_dict[k]) / len(internal_i), len(j_clade_dict[k]) / len(internal_j)]) * 50
            df.append([x, y, "Joint", s])
            if x - y > 1:
                count_outsides += len(i_clade_dict[k])
                count_outsides += len(j_clade_dict[k])

            x_var = np.std(i_clade_dict[k])
            y_var = np.std(j_clade_dict[k])
            line_coords.append([[x - (x_var / 2), x + (x_var / 2)], [y, y]])
            line_coords.append([[x, x], [y - (y_var / 2), y + (y_var / 2)]])
        else:
            if k in i_clade_dict:
                x = np.mean(i_clade_dict[k])
                y = -0.1
                s = len(i_clade_dict[k]) / len(internal_i)

                count_outsides += len(i_clade_dict[k])

                x_var = np.std(i_clade_dict[k])
                line_coords.append([[x - (x_var / 2), x + (x_var / 2)], [y, y]])
            else:
                y = np.mean(j_clade_dict[k])
                x = -0.1
                s = len(j_clade_dict[k]) / len(internal_j)

                count_outsides += len(j_clade_dict[k])

                y_var = np.std(j_clade_dict[k])
                line_coords.append([[x, x], [y - (y_var / 2), y + (y_var / 2)]])

            df.append([x, y, "Unique", s * 50])

    df = pd.DataFrame(df, columns=["X", "Y", "Type", "Size"])

    divisor = np.sum([len(v) for v in i_clade_dict.values()]) + np.sum([len(v) for v in j_clade_dict.values()])
    percentage_agreement_clades = (1 - (count_outsides / divisor)) * 100

    # todo measure the distribution of the clades ranks, if not every rank is taken then that is not good
    #  how close is the distribution to uniform, i.e. one clade per rank?
    # todo idea overlay a centroid to this plot and see what it decides the clades are!?

    if plot:
        fig, ax = plt.subplots()

        # sns.scatterplot(data=df, x="X", y="Y", s=df["Size"], ax=ax, zorder=25)
        p = sns.jointplot(data=df, x="X", y="Y", s=df["Size"], ax=ax, zorder=25,
                          marginal_kws=dict(bins=len(internal_i[0]) - 2, fill=False, discrete=True,
                                            weights=df["Size"]))
        p.ax_joint.grid()

        # adding std confidence interval lines
        for lc in line_coords:
            p.ax_joint.plot(lc[0], lc[1], color="red", linewidth=.3, zorder=1)

        # adding 1 rank deviation diagonals
        p.ax_joint.plot([0, len(internal_i[0]) - 2 - 1], [1, len(internal_i[0]) - 2], color="grey", linewidth=.5,
                        linestyle="--", zorder=20)
        p.ax_joint.plot([1, len(internal_i[0]) - 2], [0, len(internal_i[0]) - 2 - 1], color="grey", linewidth=.5,
                        linestyle="--", zorder=20)

        # setting x and y axis limits
        p.ax_joint.set_xlim(-0.75, len(internal_i[0]) - 1.5)
        p.ax_joint.set_ylim(-0.75, len(internal_i[0]) - 1.5)

        # plt.suptitle(f"{percentage_agreement_clades}")

        p.ax_joint.set_xlabel(f"{percentage_agreement_clades}")
        # p.ax_joint.set_xlabel("Rank")
        # p.ax_joint.set_ylabel("Rank")

        if file:
            plt.savefig(file, dpi=300, bbox_inches="tight")
        else:
            plt.show()
        plt.clf()
        plt.close("all")
    return percentage_agreement_clades
