import matplotlib.pyplot as plt
import scipy.stats as stat
from pingouin import multivariate_normality
import pandas as pd

class TripleNotAClade(Exception):
    pass


class TripleInNoTree(Exception):
    pass


def get_all_triples(tree):
    # Take an ete3 tree and extract all triples
    triples = []

    for node in tree.traverse('levelorder'):  # instead of recursion
        if len(node) > 2:
            # more than 2 leafs beneath the node
            c = node.children
            c0_leafs = []
            for leaf in c[0]:
                c0_leafs.append(leaf.name)
            c1_leafs = []
            for leaf in c[1]:
                c1_leafs.append(leaf.name)
            # all triples c0| c10,c11
            if len(c[0]) == 1:
                t2 = []
                for leaf in c[1].children[0]:
                    t2.append(leaf.name)
                t3 = []
                for leaf in c[1].children[1]:
                    t3.append(leaf.name)
                triples.append([c0_leafs,t2,t3])
                continue
            if len(c[1]) == 1:
                t2 = []
                for leaf in c[0].children[0]:
                    t2.append(leaf.name)
                t3 = []
                for leaf in c[0].children[1]:
                    t3.append(leaf.name)
                triples.append([c1_leafs, t2, t3])
                continue

            t2 = []
            for leaf in c[1].children[0]:
                t2.append(leaf.name)
            t3 = []
            for leaf in c[1].children[1]:
                t3.append(leaf.name)
            triples.append([c0_leafs, t2, t3])

            t2 = []
            for leaf in c[0].children[0]:
                t2.append(leaf.name)
            t3 = []
            for leaf in c[0].children[1]:
                t3.append(leaf.name)
            triples.append([c1_leafs, t2, t3])
    return triples


def get_tree_coordinates(tree, triple_leaves):

    combined_triple_set = [taxa for sublist in triple_leaves for taxa in sublist]

    # if tree.etree.get_common_ancestor(combined_triple_set).is_root():
    #    raise TripleNotASubtree

    # triple_root = tree.etree.get_common_ancestor(combined_triple_set)
    triple_root = tree.etree.get_common_ancestor(combined_triple_set)
    if len(triple_root) != len(triple_leaves[0] + triple_leaves[1] + triple_leaves[2]):
        raise TripleNotAClade

    # todo other checks
    a = triple_leaves[0]
    b = triple_leaves[1]
    c = triple_leaves[2]

    ab = triple_leaves[0] + triple_leaves[1]
    ac = triple_leaves[0] + triple_leaves[2]
    bc = triple_leaves[1] + triple_leaves[2]

    # todo switch to proper parameters x and y

    mrca_ab = tree.etree.get_common_ancestor(ab)
    mrca_ac = tree.etree.get_common_ancestor(ac)
    mrca_bc = tree.etree.get_common_ancestor(bc)

    if len(a) != 1 and len(tree.etree.get_common_ancestor(a)) != len(a):
        raise TripleNotAClade
    if len(b) != 1 and len(tree.etree.get_common_ancestor(b)) != len(b):
        raise TripleNotAClade
    if len(c) != 1 and len(tree.etree.get_common_ancestor(c)) != len(c):
        raise TripleNotAClade


    # todo currently the return is only height of the triple and the one internal branch length
    # todo triple_string should contain the actual taxa?

    # the name contains the distance of a node to the root
    if float(mrca_ab.name) > float(mrca_ac.name):
        # ab|c is the triple
        triple_string = "ab|c"
        param_y = mrca_ab.get_distance(mrca_ab.children[0])  # height of the cherry
        # triple_height = mrca_ac.get_distance(tree.etree.get_leaves_by_name(ab[0])[0])  # distance from the triple root to any leaf is the height
        param_x = mrca_ac.get_distance(mrca_ab)  # branch length root to cherry
    elif float(mrca_ac.name) > float(mrca_bc.name):
        # ac|b is the triple
        triple_string = "ac|b"
        param_y = mrca_ac.get_distance(mrca_ac.children[0])  # height of cherry
        # triple_height = mrca_ab.get_distance(tree.etree.get_leaves_by_name(ab[0])[0])  # distance from the triple root to any leaf is the height
        param_x = mrca_ab.get_distance(mrca_ac)  # branch length root to cherry
    else:
        # bc|a is the triple
        triple_string = "bc|a"
        param_y = mrca_bc.get_distance(mrca_bc.children[0])  # height of cherry
        # triple_height = mrca_ab.get_distance(tree.etree.get_leaves_by_name(ab[0])[0])  # distance from the triple root to any leaf is the height
        param_x = mrca_ab.get_distance(mrca_bc)  # branch length root to cherry

    return triple_string, param_x, param_y


def get_trees_embedding(trees, triple_leaves, tree_space):
    # dictionary containing the coordinate paris for each page
    page_coords = {"ab|c": [], "ac|b": [], "bc|a": []}
    fails = 0
    for t in trees:
        try:
            t, x, y = get_tree_coordinates(t, triple_leaves)
            if tree_space == 1:
                page_coords[t].append((x, y))
            if tree_space == 2:
                page_coords[t].append((x+y, y))
        except TripleNotAClade:
            fails += 1
            # todo count the number of trees containing the subtree
            pass
    return page_coords, fails


import numpy as np
from scipy.stats import multivariate_normal
import seaborn as sns


def plot_pages(page1, page1_label, page2, page2_label, pages_ax):
    # plotting page 1 and 2 accordingly

    _p1 = False
    _p2 = False

    x1, x2, y1, y2 = (), (), (), ()
    # page2 coordinates need to be transformed, bl * -1, because height is the shared axis, bl is the 2nd coordinate
    try:
        x2, y2 = zip(*page2)
        pages_ax.scatter([-tmp for tmp in x2], y2, label=page2_label, alpha=0.5, s=1)
        if len(x2) > 20:
            sns.kdeplot(x=[-tmp for tmp in x2], y=y2, ax=pages_ax, warn_singular=False)

        # grid_x = np.linspace(np.min(page2), np.max(page2), 10, endpoint=False)
        # grid_x, grid_y = np.mgrid[-np.max(x2):-np.min(x2):0.01, np.min(y2):np.max(y2):0.01]
        # normal = multivariate_normal([-np.mean(x2), np.mean(y2)], [np.cov(x2), np.cov(y2)])
        # pos = np.dstack((grid_x, grid_y))
        # pages_ax.contour(grid_x, grid_y, normal.pdf(pos))
        _p2 = True
    except ValueError:
        pass

    try:
        x1, y1 = zip(*page1)
        pages_ax.scatter(x1, y1, label=page1_label, alpha=0.5, s=1)
        if len(x1) > 20:
            sns.kdeplot(x=x1, y=y1, ax=pages_ax, warn_singular=False)

        # grid_x, grid_y = np.mgrid[np.min(x1):np.max(x1):0.01, np.min(y1):np.max(y1):0.01]
        # normal = multivariate_normal([np.mean(x1), np.mean(y1)], [np.cov(x1), np.cov(y1)])
        # pos = np.dstack((grid_x, grid_y))
        # pages_ax.contour(grid_x, grid_y, normal.pdf(pos))
        _p1 = True
    except ValueError:
        pass

    t1_pingouin, t2_pingouin = 0, 0
    if len(x1) > 19:
        # t_scipy = stat.normaltest((x1, y1), axis=1)
        t1_pingouin = multivariate_normality(pd.DataFrame({"x": x1, "y": y1})).pval

    if len(x2) > 19:
        # t_scipy = stat.normaltest((x1, y1), axis=1)
        t2_pingouin = multivariate_normality(pd.DataFrame({"x": x2, "y": y2})).pval

    if _p1 or _p2:
        pages_ax.set_title(f"{t2_pingouin}, {t1_pingouin}")
        # pages_ax.setylim(0, max([max(y1), max(y2)])+.1)
        pages_ax.axvline(x=0, color="black")
        pages_ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.05),
                        ncol=3, fancybox=True, shadow=True)

    return 0


def distribution_classifier(trees, triples, tree_space=1, page_coords=False, fails=False):

    if not page_coords or fails:
        page_coords, fails = get_trees_embedding(trees, triples, tree_space=tree_space)

    # if not (1 - (fails / len(trees))) > .5:
    #     return None

    ab = len(page_coords["ab|c"])
    ac = len(page_coords["ac|b"])
    bc = len(page_coords["bc|a"])
    total_points = ab + ac + bc
    if total_points == 0:
        raise TripleInNoTree()
    values = sorted([ab/total_points, ac/total_points, bc/total_points], reverse=True)

    if values[0] >= .95:
        return "A", (1- (fails/len(trees))), (values[0], values[1])
    if values[0] + values[1] >= .95:
        if values[1] >= (values[0]/2):
            return "AB", (1- (fails/len(trees))), (values[0], values[1])
        return "Ab", (1- (fails/len(trees))), (values[0], values[1])
    if values[1] >= (values[0] / 2):
        if values[2] >= (values[0] / 2):
            return "ABC", (1- (fails/len(trees))), (values[0], values[1])
        return "ABc", (1- (fails/len(trees))), (values[0], values[1])
    return "Abc", (1- (fails/len(trees))), (values[0], values[1])


def plot_book(trees, triples, tree_space=1, add_class=False):

    fig, axs = plt.subplots(3, sharex=True, figsize=(12, 8))

    page_coords, fails = get_trees_embedding(trees, triples, tree_space=tree_space)
    try:
        dc, posterior_support, _ = distribution_classifier(trees, triples, tree_space=tree_space, page_coords=page_coords, fails=fails)
    except TripleInNoTree:
        dc = ""
        pass

    plt.suptitle(f"%trees = {(1 - (fails/len(trees)))*100}\na={triples[0]}, b={triples[1]}, c={triples[2]}\n{dc if add_class else ''}")

    plot_pages(page1=page_coords["ab|c"], page2=page_coords["bc|a"], page1_label="ab|c", page2_label="bc|a", pages_ax=axs[0])
    plot_pages(page1=page_coords["ab|c"], page2=page_coords["ac|b"], page1_label="ab|c", page2_label="ac|b", pages_ax=axs[1])
    plot_pages(page1=page_coords["bc|a"], page2=page_coords["ac|b"], page1_label="bc|a", page2_label="ac|b", pages_ax=axs[2])

    return fig
    # plt.show()
