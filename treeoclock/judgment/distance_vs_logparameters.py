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


def plot_log_neighbours(mchain, threshold):
    # todo a proper test case for this is necessary if i develop more about this part

    log_key = "posterior"

    distances = mchain.pwd_matrix()
    log_values = mchain.log_data[log_key]

    if distances.shape[0] < log_values.shape[0] and (mchain.tree_sampling_interval % mchain.log_sampling_interval == 0):
        sampling_diff = int(mchain.tree_sampling_interval / mchain.log_sampling_interval)
        log_values = log_values[::2]  # deleting logs that do not have a tree
    elif distances.shape[0] > log_values.shape[0]:
        raise ValueError("Problem with tree sample size!")

    log_values = np.array(log_values)

    # todo delete multiples of the same tree, for this delete from matrix and from log values list

    max_diff = np.max(log_values) - np.min(log_values)

    seen = {}
    df = []
    for i in range(distances.shape[0]):
        l1 = np.where((distances[i, :] < threshold) & (distances[i, :] > 0))[0]
        # l1 = np.append(l1, np.where(distances[:, i] == 1)[0])
        if np.all(distances[:(i-1), i]):
            for x in l1:
                if np.all(distances[:(x-1), x]):
                    # df.append([log_values[i], log_values[x]])
                    df.append([abs(log_values[i] - log_values[x])/max_diff])
    df = pd.DataFrame(df, columns=["Value"])
    sns.violinplot(data=df, x="Value", inner="stick", cut=0)
    # df = pd.DataFrame(df, columns=["x", "y"])
    # sns.scatterplot(data=df, x="x", y="y")

    # ax = plt.axis()  # Getting current axis limits
    # plt.axline([0, 0], [1, 1], color='red')  # Plot f(x)=x
    # plt.axline([0, max_diff], [1, 1+max_diff], color="orange")
    # plt.axline([0, -max_diff], [1, 1-max_diff], color="orange")
    # plt.axis(ax)  # Resetting the axis limits

    # todo make the plot look nice and add labels and label sizes for the paper!

    plt.show()

    return 0