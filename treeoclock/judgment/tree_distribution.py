from treeoclock.summary.centroid import Centroid
import matplotlib.pyplot as plt
import seaborn as sns
from sympy import Sum, Product
from sympy.abc import x, i, m
from treeoclock.judgment.conditional_clade_distribution import get_tree_probability, get_maps, add_centroid
import pandas as pd
import numpy as np
import scipy.stats as st


def get_coefficients(n):
    # i, m = symbols('i m', integer=True)
    h = n-2
    c = Product(Sum(i*x**(i-1), (i, 1, m)), (m, 1, h+1)).doit().as_poly().coeffs()
    c.reverse()
    return [int(co) for co in c]


def _get_approx_orbits(n):
    # from wolframclient.language import wlexpr
    # from wolframclient.evaluation import WolframLanguageSession
    # session = WolframLanguageSession()
    # # the last loop for n is number of taxa -2
    # expression = f'Table[CoefficientList[Product[Sum[D[x^i, x], {{i, 0, m + 1}}], {{m, 0, n}}], x], {{n, {n - 2},{n - 2}}}]'
    # r = session.evaluate(wlexpr(expression))
    # coeff = list(r[0])
    # session.terminate()
    # return coeff

    # first calculate the sum
    sumands = []
    for i in range(1,n):
        sumands.append((i, i-1))


    # calculate the product
    product = {0: 1}  # initialized with the trees at distance 0
    for k in range(1,n):
        # tmp_product = product.copy()
        high_prev = max(product.keys())
        tmp_product = {i: 1 for i in range(high_prev+k)}
        # print(tmp_product.keys(), tmp_product.items())
        for coeff_s, exp_s in sumands[0:k]:
            for exp, coeff in product.items():
                if tmp_product[exp + exp_s] == 1:
                    tmp_product[exp + exp_s] = coeff * coeff_s
                else:
                    tmp_product[exp + exp_s] += (coeff * coeff_s)
        product = tmp_product.copy()
    return list(product.values())


def plot_tree_density_distribution(Mchain, centroid="calc", given_x=-1, ix_chain=0):
    # todo the MChain object has a centroid object in it, use that one at some point
    if centroid == "calc":
        mycen = Centroid(variation="inc_sub", n_cores=24)
        centroid, sos = mycen.compute_centroid(Mchain[ix_chain].trees)

    cen_distances = [t.fp_distance(centroid) for t in Mchain[ix_chain].trees]

    n = len(centroid)
    if given_x == -1:
        # given_x = max(cen_distances)
        given_x = int(((n-1)*(n-2))/2)

        # orbits = _get_approx_orbits(len(centroid))
    orbits = get_coefficients(n)  # total number of trees at distance
    points = []  # number of trees at distanes to the centroid
    points_2 = []  # fraction of orbit at distance to centroid that has been sampled
    points_3 = []  # samples at distance from centroid
    # for i in range(len(orbits)-1):
    for i in range(given_x):
        points.append(cen_distances.count(i)/orbits[i])
        points_2.append(sum(d < i for d in cen_distances)/sum(orbits[0:i+1]))
        points_3.append(cen_distances.count(i))
    fig, ax = plt.subplots(2,2, sharex=True)

    # sns.lineplot(x = range(len(orbits)-1), y = points, ax=ax[0])
    sns.lineplot(x = range(given_x), y = points, ax=ax[0, 0])

    # sns.lineplot(x=range(len(orbits) - 1), y=points_2, ax=ax[1])
    sns.lineplot(x=range(given_x), y=points_2, ax=ax[1, 0])


    sns.lineplot(x=range(given_x), y=points_3, ax=ax[0, 1])

    sns.lineplot(x=range(len(orbits)), y=orbits, ax=ax[1, 1])

    # print(orbits[0:min(max(cen_distances), given_x)+1])

    ax[1, 0].set_xlabel("Distance from centroid")
    ax[0, 0].set_ylabel("Ts at d/\n orbit at d")
    ax[1, 0].set_ylabel("sum of Ts upto d/\n sum orbits upto d")

    ax[0, 0].set_title(f"{n} taxa, max distance = {int(((n-1)*(n-2))/2)}")
    fig.show()
    return 0


def get_sample_treespace_coverage(Mchain, centroid="calc", ix_chain=0):

    # todo this is probably way to inefficient, but it'll do for now

    if centroid == "calc":
        mycen = Centroid(variation="inc_sub", n_cores=24)
        centroid, _ = mycen.compute_centroid(Mchain[ix_chain].trees)

    cen_distances = [t.fp_distance(centroid) for t in Mchain[ix_chain].trees]
    points = []
    n = len(centroid)
    for i in range(int(((n-1)*(n-2))/2)):
        points.append(cen_distances.count(i))
    orbits = get_coefficients(n)  # total number of trees at distance
    
    total_trees = len(Mchain[ix_chain].trees)
    coverages = [sum(points[0:i+1]) / total_trees for i in range(len(points))]

    index_95 = -1
    for ix, el in enumerate(coverages):
        if el >= 0.95:
            index_95 = ix
            break
    return sum(points[0:index_95+1])/sum(orbits[0:index_95+1])


def plot_CCD_vs_centroid_distance(Mchain, ix_chain = 0, centroid = "calc"):

    if centroid == "calc":
        mycen = Centroid(variation="inc_sub", n_cores=24)
        centroid, _ = mycen.compute_centroid(Mchain[ix_chain].trees)

    m1, m2, uniques = get_maps(Mchain[ix_chain].trees)

    # todo work with the new uniques dictionary to calculate the correct plot


    cen_distances = [t.fp_distance(centroid) for t in Mchain[ix_chain].trees]
    points = []
    n = len(centroid)
    for i in range(int(((n - 1) * (n - 2)) / 2)):
        points.append(cen_distances.count(i))

    # total_trees = len(Mchain[ix_chain].trees)
    # coverages = [sum(points[0:i + 1]) / total_trees for i in range(len(points))]
    # finding 95 percent radius orbit, i.e. at this distance 95% of samples are closer to the tree
    # index_95 = -1
    # for ix, el in enumerate(coverages):
    #     if el >= 0.95:
    #         index_95 = ix
    #         break
    index_95 = np.max(cen_distances)

    data = []

    # todo different?
    cen_probability = get_tree_probability(centroid, m1, m2)

    # if cen_probability == 0:
    #     # todo this is still an open problem i suppose
    #     m1, m2 = add_centroid(centroid, m1, m2)
    # cen_probability = get_tree_probability(centroid, m1, m2)

    for t in uniques.keys():
    # for t in range(len(Mchain[ix_chain].trees)):
        cur_probability = get_tree_probability(Mchain[ix_chain].trees[t], m1, m2)
        mean_cen_distance = np.mean([cen_distances[k] for k in uniques[t]] + [cen_distances[t]])
        if mean_cen_distance <= index_95:  # all trees that are within the 95% radius of centroid
            data.append([mean_cen_distance, np.log(cur_probability), len(uniques[t]) == 0])
    data = pd.DataFrame(data, columns = ["Dist", "Prob", "Shift"])

    # todo fit linear regression to the data with logscale to the probability value!
    #  check the sos ll correlation for how to get the correlation value of these
    #  This is going to just be correlation for now!
    # todo

    # pearson correlation changes with log
    pr = st.pearsonr(data["Prob"], data["Dist"])

    data_f = data.query("Shift == False")
    data_t = data.query("Shift == True")

    pr_f = "Not enough samples"
    pr_t = "Not enough samples"
    if data_f.shape[0] > 5:
        pr_f = st.pearsonr(data_f["Prob"], data_f["Dist"])
    if data_t.shape[0] > 5:
        pr_t = st.pearsonr(data_t["Prob"], data_t["Dist"])

    # spearman is invariant to monotone transformation becausse based on ranks
    # sr = st.spearmanr(np.log(data["Prob"]), data["Dist"])


    # bp = sns.boxplot(data=data, x="Dist", y="Prob", order=sorted(set(data["Dist"].values)))
    # bp = sns.scatterplot(data=data, x="Dist", y="Prob")
    bp = sns.lmplot(x="Dist", y="Prob", data=data, hue="Shift", legend=False)
    # bp = sns.regplot(y="Dist", x="Prob", data=data, logx=True, line_kws={"color": "red"})
    
    # xmin, xmax = bp.get_ylim()
    # bp.vlines(x=cen_probability, ymin=xmin, ymax=xmax, ls="--", lw=4, colors="tab:orange")
    # plt.axhline(y = cen_probability, color="purple")

    plt.legend(loc="lower left")

    plt.ylabel("CCD (Probability), Log Scale")
    plt.suptitle(f"Correlation: {pr[0]}\nMultiple(false): {pr_f[0]}\nSingle(true): {pr_t[0]}")
    plt.xlabel("Distance to centroid")

    # plt.xscale('log')

    plt.xticks(rotation=270)
    plt.tight_layout()
    # plt.show()
    # return 0

    plt.savefig(
        fname=f"{Mchain.working_dir}/plots/{Mchain.name}_{ix_chain}_cen_dist_CCD.png",
        format="png", bbox_inches="tight", dpi=800)
    plt.clf()
    plt.close("all")
    return 0
