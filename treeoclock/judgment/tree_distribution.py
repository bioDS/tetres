from treeoclock.summary.centroid import Centroid
import matplotlib.pyplot as plt
import seaborn as sns
from sympy import Sum, Product
from sympy.abc import x, i, m
from treeoclock.judgment.conditional_clade_distribution import get_tree_probability, get_maps

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


def get_fraction_95area():
    # todo get the 95% interval of the trees and then look at the fraction with all trees vs. up to that orbit size

    return 1


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
        centroid, sos = mycen.compute_centroid(Mchain[ix_chain].trees)

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


def plot_CCD_vs_centroid_distance(Mchain, ix_chain = 0):

    mycen = Centroid(variation="inc_sub", n_cores=24)
    centroid, sos = mycen.compute_centroid(Mchain[ix_chain].trees)
    m1, m2, uniques = get_maps(Mchain[ix_chain].trees)

    cen_distances = [Mchain[ix_chain].trees[t].fp_distance(centroid) for t in uniques]
    tree_probs = [get_tree_probability(Mchain[ix_chain].trees[t], m1, m2) for t in uniques]

    # cen_distances = [t.fp_distance(centroid) for t in Mchain[0].trees]
    # tree_probs = [get_tree_probability(t, m1, m2) for t in Mchain[0].trees]

    sns.boxplot(x=cen_distances, y=tree_probs)
    plt.ylabel("CCD (Probability)")
    plt.suptitle("Comparing Larget CCD probaility with distance to centroid tree")
    plt.xlabel("Distance to centroid")

    plt.yscale('log')

    plt.savefig(
        fname=f"{Mchain.working_dir}/plots/{Mchain.name}_{ix_chain}_cen_dist_CCD.png",
        format="png", bbox_inches="tight", dpi=800)
    plt.clf()
    plt.close("all")

    return 1
