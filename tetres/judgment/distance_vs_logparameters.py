import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# pip install git+https://github.com/widmi/multiprocess-shared-numpy-arrays
from share_array.share_array import get_shared_array, make_shared_array
from multiprocessing import Pool
import itertools
from rpy2.robjects.packages import importr


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


# todo add to the pairwise distance matrix code with saving and all that instead of here ....
def _kc(t1, t2, i, j, kc_param = 0):
    matrix = get_shared_array('distances')
    ape = importr("ape")
    # TreeDist = importr("TreeDist")
    treespace = importr("treespace")
    # distance = TreeDist.KendallColijn(ape.read_tree(text = t1), ape.read_tree(text = t2))[0]
    distance = treespace.treeDist(ape.read_tree(text = t1), ape.read_tree(text = t2), kc_param)[0]
    matrix[i, j] = distance
    return 0


def get_kc_matrix(trees, kc_param):
    n = len(trees)
    distances = np.zeros((n, n))
    make_shared_array(distances, name='distances')  # create shared memory array from numpy array
    # shared_array = get_shared_array('distances')  # get shared memory array as numpy array
    with Pool(60) as p:
        p.starmap(_kc, [(trees[i].get_newick(), trees[j].get_newick(), i, j, kc_param) for i, j in itertools.combinations(range(n), 2)])
    distances = get_shared_array('distances')
    try:
        del globals()['distances']
    except KeyError:
        pass
    return distances


def plot_log_neighbours(mchain, dist_type = "rnni", kc_param=0):

    if dist_type not in ["rnni", "rf", "kc"]:
        raise ValueError(f"Distance {dist_type} not supported!")

    df = []
    for log_key in ["posterior", "likelihood"]:
        if dist_type == "rf":
            distances = mchain.pwd_matrix(rf=True)
        elif dist_type == "rnni":
            distances = mchain.pwd_matrix(rf=False)
        elif dist_type == "kc":
            distances = get_kc_matrix(mchain.trees, kc_param=kc_param)
            # todo

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
            l1 = range(i+1, distances.shape[0])
            # l1 = np.where((distances[i, :] < int(0.25 * max_d)) | (distances[i, :] > int(0.75 * max_d)))[0]
            for x in l1:
                if np.all(distances[:(i - 1), i]):
                    if x != i:
                        if np.all(distances[:(x - 1), x]):
                            if round(distances[i, x]/max_d, 7) != 0:
                                df.append([abs(log_values[i] - log_values[x])/max_diff, log_key, round(distances[i, x]/max_d, 7)])


    df = pd.DataFrame(df, columns=["Value", "Parameter", "Offset"])
    print("Plotting")
    sns.lineplot(data=df, y="Value", x="Offset", hue="Parameter")
    print("Finished")

    plt.ylabel("Diff. of log parameter\n (relative)", fontsize=16)
    plt.xlabel(f"{dist_type} distance (relative)", fontsize=16)

    plt.tick_params(labelsize=14)
    plt.tight_layout()

    plt.savefig(f"{mchain.working_dir}/plots/smoothness_plot_{dist_type}_{kc_param if dist_type == 'kc' else ''}.eps", format="eps", dpi=400, bbox_inches="tight")
    # plt.show()
    plt.clf()
    plt.close()
    return 0
