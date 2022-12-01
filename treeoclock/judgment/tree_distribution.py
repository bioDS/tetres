from treeoclock.summary.centroid import Centroid
import matplotlib.pyplot as plt
import seaborn as sns
from sympy import Sum, Product
from sympy.abc import x, i, m


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


def plot_tree_density_distribution(Mchain, centroid="calc", given_x=-1):
    # todo the MChain object has a centroid object in it, use that one at some point
    if centroid == "calc":
        mycen = Centroid(variation="inc_sub", n_cores=24)
        centroid, sos = mycen.compute_centroid(Mchain[0].trees)

    cen_distances = [t.fp_distance(centroid) for t in Mchain[0].trees]

    if given_x == -1:
        given_x = max(cen_distances)

    # orbits = _get_approx_orbits(len(centroid))
    orbits = get_coefficients(len(centroid))
    n = len(centroid)
    points = []
    points_2 = []
    points_3 = []
    # for i in range(len(orbits)-1):
    for i in range(min(max(cen_distances), given_x)):
        points.append(cen_distances.count(i)/orbits[i])
        points_2.append(sum(d < i for d in cen_distances)/sum(orbits[0:i+1]))
        points_3.append(cen_distances.count(i))
    fig, ax = plt.subplots(2,2, sharex=True)

    # sns.lineplot(x = range(len(orbits)-1), y = points, ax=ax[0])
    sns.lineplot(x = range(min(max(cen_distances), given_x)), y = points, ax=ax[0, 0])

    # sns.lineplot(x=range(len(orbits) - 1), y=points_2, ax=ax[1])
    sns.lineplot(x=range(min(max(cen_distances), given_x)), y=points_2, ax=ax[1, 0])

    
    sns.lineplot(x=range(min(max(cen_distances), given_x)), y=points_3, ax=ax[0, 1])
    # sns.lineplot(x=range(min(max(cen_distances), given_x)), y=orbits[0:min(max(cen_distances), given_x)+1], ax=ax[0, 1])


    ax[1, 0].set_xlabel("Distance from centroid")
    ax[0, 0].set_ylabel("Ts at d/\n orbit at d")
    ax[1, 0].set_ylabel("sum of Ts upto d/\n sum orbits upto d")

    ax[0, 0].set_title(f"{n} taxa, max distance = {int(((n-1)*(n-2))/2)}")

    fig.show()

    return 0