import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# todo make this a more elaborate thing in the future?


def plot_distance_posteroir_comparison(mchain):
    log_key = "posterior"

    distances = mchain.pwd_matrix()
    log_values = mchain.log_data[log_key]

    if distances.shape[0] < log_values.shape[0] and (mchain.tree_sampling_interval % mchain.log_sampling_interval == 0):
        sampling_diff = int(mchain.tree_sampling_interval / mchain.log_sampling_interval)
        log_values = log_values[::2]  # deleting logs that do not have a tree
    elif distances.shape[0] > log_values.shape[0]:
        raise ValueError("Problem with tree sample size!")

    log_values = np.array(log_values)

    plot_data = []

    for i in range(1, distances.shape[0]-1):
        for j in range(i+1, distances.shape[0]):
            cur_log_diff = abs(log_values[i]-log_values[j])
            plot_data.append([log_key, distances[i, j], cur_log_diff])

    plot_data = pd.DataFrame(plot_data, columns=["Log", "RNNI distance", f"{log_key} difference"])
    
    # sns.scatterplot(data=plot_data, x="RNNI distance", y=f"{log_key} difference")
    sns.scatterplot(data = plot_data.nlargest(n=100, columns=f"{log_key} difference"),
                    x="RNNI distance", y=f"{log_key} difference", color="red")
    sns.scatterplot(data=plot_data.nsmallest(n=100, columns=f"{log_key} difference"),
                    x="RNNI distance", y=f"{log_key} difference", color="blue")
    plt.show()
    plt.clf()


def plot_sos_log_comparison(mchain):
    log_key = "posterior"

    distances = mchain.pwd_matrix()
    log_values = mchain.log_data[log_key]

    if distances.shape[0] < log_values.shape[0] and (mchain.tree_sampling_interval % mchain.log_sampling_interval == 0):
        sampling_diff = int(mchain.tree_sampling_interval / mchain.log_sampling_interval)
        log_values = log_values[::2]  # deleting logs that do not have a tree
    elif distances.shape[0] > log_values.shape[0]:
        raise ValueError("Problem with tree sample size!")

    log_values = np.array(log_values)

    plot_data = []
    for i in range(distances.shape[0]):
        cur_log = log_values[i]
        plot_data.append([log_key, np.sum(distances[i, :]**2)+np.sum(distances[:, i]**2), cur_log])

    plot_data = pd.DataFrame(plot_data, columns=["Log", "SoS", f"{log_key}"])
    plot_data.drop(plot_data.nsmallest(n=10, columns=f"{log_key}").index, inplace=True)

    # sns.scatterplot(data=plot_data, x="SoS", y=f"{log_key}")
    sns.scatterplot(data=plot_data.nlargest(n=100, columns=f"{log_key}"),
                    x="SoS", y=f"{log_key}", color="red")
    sns.scatterplot(data=plot_data.nsmallest(n=100, columns=f"{log_key}"),
                    x="SoS", y=f"{log_key}", color="blue")
    plt.show()
    plt.clf()


def plot_loglik_along_path(mchain):
    log_key = "likelihood"

    distances = mchain.pwd_matrix()
    log_values = mchain.log_data[log_key]

    if distances.shape[0] < log_values.shape[0] and (mchain.tree_sampling_interval % mchain.log_sampling_interval == 0):
        sampling_diff = int(mchain.tree_sampling_interval / mchain.log_sampling_interval)
        log_values = log_values[::2]  # deleting logs that do not have a tree
    elif distances.shape[0] > log_values.shape[0]:
        raise ValueError("Problem with tree sample size!")

    log_values = np.array(log_values)

    x = 0
    df = []
    for i in range(distances.shape[0]):
        df.append([x, log_values[i]])
        if i < distances.shape[0]-1:
            x = distances[0, i]

    df = pd.DataFrame(df, columns=["Sample", f"{log_key}"])
    sns.lineplot(data=df, x="Sample", y=f"{log_key}", marker="o")
    plt.xlabel("RNNI dist to tree at x = 0")
    plt.show()
    plt.clf()
    return 0


def plot_log_neighbours(mchain):

    df = []
    for log_key in ["posterior", "likelihood"]:
        distances = mchain.pwd_matrix()
        log_values = mchain.log_data[log_key]

        if distances.shape[0] < log_values.shape[0] and (mchain.tree_sampling_interval % mchain.log_sampling_interval == 0):
            sampling_diff = int(mchain.tree_sampling_interval / mchain.log_sampling_interval)
            log_values = log_values[::sampling_diff]  # deleting logs that do not have a tree
        elif distances.shape[0] > log_values.shape[0]:
            raise ValueError("Problem with tree sample size!")

        log_values = np.array(log_values)

        max_diff = abs(np.max(log_values) - np.min(log_values))
        max_d = np.max(distances)
        for i in range(distances.shape[0]):
            # low distances
            l1 = np.where((distances[i, :] < int(0.1*max_d)) & (distances[i, :] > 0))[0]
            # l1 = np.append(l1, np.where(distances[:, i] == 1)[0])
            if np.all(distances[:(i-1), i]):
                for x in l1:
                    if np.all(distances[:(x-1), x]):
                        # df.append([log_values[i], log_values[x]])
                        df.append([abs(log_values[i] - log_values[x])/max_diff, log_key, f"0-{int(0.1*max_d)}"])
            # high distances
            l1 = np.where((distances[i, :] < max_d) & (distances[i, :] > int(0.9*max_d)))[0]
            # l1 = np.append(l1, np.where(distances[:, i] == 1)[0])
            if np.all(distances[:(i-1), i]):
                for x in l1:
                    if np.all(distances[:(x-1), x]):
                        # df.append([log_values[i], log_values[x]])
                        df.append([abs(log_values[i] - log_values[x])/max_diff, log_key, f"{int(0.9*max_d)}-{max_d}"])
    df = pd.DataFrame(df, columns=["Value", "Parameter", "Offset"])
    # sns.violinplot(data=df, x="Value", y="Parameter", inner="stick", cut=0)
    sns.boxplot(data=df, y="Value", x="Parameter", hue="Offset")
    # plt.xlim(0, 1)
    plt.ylabel("Difference of normalized log value")
    low_samples = len(df[(df["Offset"] == f"0-{int(0.1*max_d)}") & (df["Parameter"] == "likelihood")])
    high_samples = len(df[(df["Offset"] == f"{int(0.9*max_d)}-{max_d}") & (df["Parameter"] == "likelihood")])
    # plt.suptitle(f"Comparing log parameters for {samples} trees with distance < {threshold}")
    plt.suptitle(f"Comparing log parameters\n10% lowest distances ({low_samples}) vs. 10% highest distances ({high_samples})")

    plt.savefig(f"{mchain.working_dir}/plots/smoothness_plot.png", dpi=400, bbox_inches="tight")
    # plt.show()
    plt.clf()
    plt.close()
    return 0