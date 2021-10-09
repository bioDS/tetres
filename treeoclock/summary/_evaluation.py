from treeoclock.summary.centroid import Centroid
from treeoclock.trees.time_trees import TimeTreeSet

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import timeit

# Code used for timing and comparing different centroid variations


def eval_plot(data, name):
    data = pd.DataFrame(data, columns=["Variation", "Type", "Value"])

    fig, axs = plt.subplots(2, 1, sharex='col')
    fig.tight_layout()

    sns.swarmplot(ax=axs[0], data=data[(data["Type"] == "sos")], x="Variation", y="Value")
    sns.swarmplot(ax=axs[1], data=data[(data["Type"] == "time")], x="Variation", y="Value")

    axs[0].set_ylabel('SOS value')
    axs[1].set_ylabel('Time in Seconds')

    plt.subplots_adjust(hspace=.0)  # removing vertical space between plots
    plt.xlabel(f'Centroid Variants')

    plt.savefig(f'{name}_eval.png', dpi=200, bbox_inches='tight')
    plt.close()
    # plt.show()
    

def eval_func(name):
    greedy = Centroid(start="FM", variation="greedy")
    separate = Centroid(start="FM", variation="separate")
    onlyone = Centroid(start="FM", variation="onlyone")

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{name}/{name}.trees')

    data = []

    for _ in range(3):
        start_time = timeit.default_timer()
        cen, sos = greedy.compute_centroid(myts)
        end_time = timeit.default_timer() - start_time
        data.append(["Greedy", "time", end_time])
        data.append(["Greedy", "sos", sos])

        start_time = timeit.default_timer()
        cen, sos = separate.compute_centroid(myts)
        end_time = timeit.default_timer() - start_time
        data.append(["Separate", "time", end_time])
        data.append(["Separate", "sos", sos])

        start_time = timeit.default_timer()
        cen, sos = onlyone.compute_centroid(myts)
        end_time = timeit.default_timer() - start_time
        data.append(["Onlyone", "time", end_time])
        data.append(["Onlyone", "sos", sos])

    eval_plot(data, name=name)


if __name__ == '__main__':
    # for d_name in ['Dengue', 'RSV2', '75_800_005_0', '100_800_005_1']:
    #     eval_func(d_name)
    eval_func('Dengue')
