import matplotlib.pyplot as plt


class TripleNotASubtree(Exception):
    pass


def get_tree_coordinates(tree, triple_leaves):

    combined_triple_set = [taxa for sublist in triple_leaves for taxa in sublist]

    if tree.etree.get_common_ancestor(combined_triple_set).is_root():
        raise TripleNotASubtree

    # triple_root = tree.etree.get_common_ancestor(combined_triple_set)

    ab = triple_leaves[0] + triple_leaves[1]
    ac = triple_leaves[0] + triple_leaves[2]
    bc = triple_leaves[1] + triple_leaves[2]

    mrca_ab = tree.etree.get_common_ancestor(ab)
    mrca_ac = tree.etree.get_common_ancestor(ac)
    mrca_bc = tree.etree.get_common_ancestor(bc)

    # todo currently the return is only height of the triple and the one internal branch length
    # todo triple_string should contain the actual taxa?

    # the name contains the distance of a node to the root
    if float(mrca_ab.name) > float(mrca_ac.name):
        # ab|c is the triple
        triple_string = "ab|c"
        triple_height = mrca_ac.get_distance(tree.etree.get_leaves_by_name(ab[0])[0])  # distance from the triple root to any leaf is the height
        triple_branch_length = mrca_ac.get_distance(mrca_ab)
    elif float(mrca_ac.name) > float(mrca_bc.name):
        # ac|b is the triple
        triple_string = "ac|b"
        triple_height = mrca_ab.get_distance(tree.etree.get_leaves_by_name(ab[0])[0])  # distance from the triple root to any leaf is the height
        triple_branch_length = mrca_ab.get_distance(mrca_ac)
    else:
        # bc|a is the triple
        triple_string = "bc|a"
        triple_height = mrca_ab.get_distance(tree.etree.get_leaves_by_name(ab[0])[0])  # distance from the triple root to any leaf is the height
        triple_branch_length = mrca_ab.get_distance(mrca_bc)

    return triple_string, triple_height, triple_branch_length


def get_trees_embedding(trees, triple_leaves):
    # dictionary containing the coordinate paris for each page
    page_coords = {"ab|c": [], "ac|b": [], "bc|a": []}
    for t in trees:
        try:
            t, h, bl = get_tree_coordinates(t, triple_leaves)
            page_coords[t].append((h, bl))
        except TripleNotASubtree:
            pass
    return page_coords


def plot_pages(page1, page1_label, page2, page2_label):
    # plotting page 1 and 2 accordingly
    fig, ax = plt.subplots(figsize=(12, 6))

    x1, y1 = zip(*page1)
    ax.scatter(x1, y1, label=page1_label)

    # page2 coordinates need to be transformed, bl * -1, because height is the shared axis, bl is the 2nd coordinate
    x2, y2 = zip(*page2)
    ax.scatter([-tmp for tmp in x2], y2, label=page2_label)

    plt.ylim(0, max([max(y1), max(y2)])+.1)
    plt.axvline(x=0, color="black")
    plt.legend()


    plt.show()



