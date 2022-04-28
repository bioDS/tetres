import seaborn as sns
import matplotlib.pyplot as plt


def _ess_trace_plot(data_frame):
    # todo the plot should contain which ess method was used
    #  this should also be possible to compute with mutliple methods then creating some matrix plot

    # todo the x axis should also show the actual sample number and not just
    #  the number of it in the logfile
    sns.lineplot(data=data_frame, x="Upper_i", y="Ess_value", hue="Ess_key")

    plt.axhline(y=200, color="red", linestyle='-')
    plt.axhline(y=500, color="red", linestyle='-')
    plt.axhline(y=625, color="red", linestyle='-')

    plt.show()
    plt.clf()
    return 0


def _log_trace_plot(data_frame):
    # To plot traces of values, should maybe only be one specific key

    sns.lineplot(data=data_frame, x="sample_i", y="Value")

    plt.show()
    plt.clf()
    return 0