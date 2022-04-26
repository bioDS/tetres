import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def _ess_trace_plot(data_frame):
    sns.lineplot(data=data_frame, x="Upper_i", y="Ess_value", hue="Ess_key")
    plt.show()
    plt.clf()
