from tetres.trees.time_trees import TimeTree, TimeTreeSet
from tetres.summary.compute_sos import compute_sos_mt, compute_hop_sos_mt
from tetres.summary._constants import SELECT_LIST
from tetres.trees._ctrees import TREE_LIST, TREE, PAIR

import ctypes, os
from ctypes import CDLL, POINTER

from random import choice


class NoBetterNeighbourFound(Exception):
    """Raised when no neighbour has a better SOS value"""
    pass


def intelligent_neighbourhood(t: TimeTree, clades):
    neighbours_int = []

    # sample all NNI neighbours that do contain all clades
    nn = t.nni_neighbours()
    for i in nn:
        cur_clades = i.get_clades()
        if all(x in cur_clades for x in clades):
            neighbours_int.append(i)

    # only compute rank neighbours if t contains all clades
    # Todo think about this, does it actually make sense?
    tclades = t.get_clades()
    if all(x in tclades for x in clades):
        neighbours_int.extend(t.rank_neighbours())

    return neighbours_int


def intelligent_neighbourhood_separate(t: TimeTree, clades, rank):
    neighbours_int = []

    if not rank:
        # sample all NNI neighbours that do contain all clades
        nn = t.nni_neighbours()
        for i in nn:
            cur_clades = i.get_clades()
            if all(x in cur_clades for x in clades):
                neighbours_int.append(i)

    else:
        # only compute rank neighbours if t contains all clades
        tclades = t.get_clades()
        if all(x in tclades for x in clades):
            neighbours_int.extend(t.rank_neighbours())

    return neighbours_int


def search_neighbourhood_greedy(t: TimeTree, trees: TimeTreeSet, t_value: int, n_cores=None, select='first'):
    if select not in SELECT_LIST:
        raise ValueError(f"The 'select' parameter should be"
                         f" in {SELECT_LIST} but {select} was given.")

    trees.get_common_clades()

    neighbourhood = intelligent_neighbourhood(t, trees.common_clades)

    better_neighbours = []
    cur_best = t_value
    for n in neighbourhood:
        sos = compute_sos_mt(t=n, trees=trees, n_cores=n_cores)
        if sos < cur_best:
            cur_best = sos
            better_neighbours = [n]
        elif sos == cur_best and sos < t_value:
            better_neighbours.append(n)

    if len(better_neighbours) == 0:
        # No better neighbour was found
        raise NoBetterNeighbourFound
    elif len(better_neighbours) == 1:
        # Only one neighbour with smallest value found
        return better_neighbours[0], cur_best

    # print("Multiple choices!")
    # multiple neighbours with same SOS value found, choosing by the given selection method
    if select == 'first':
        return better_neighbours[0], cur_best
    elif select == 'last':
        return better_neighbours[-1], cur_best
    elif select == 'random':
        return choice(better_neighbours), cur_best


def search_neighbourhood_greedy_omp(t: TimeTree, trees: TimeTreeSet, t_value: int, n_cores=-1, select='first'):
    if select not in SELECT_LIST:
        raise ValueError(f"The 'select' parameter should be"
                         f" in {SELECT_LIST} but {select} was given.")
    if n_cores is None:
        n_cores = -1

    trees.get_common_clades()

    neighbourhood = intelligent_neighbourhood(t, trees.common_clades)

    lib = CDLL(f"{os.path.dirname(os.path.dirname(os.path.realpath(__file__)))}/trees/findpath.so")

    lib.best_sos_neighbour.argtypes = [POINTER(TREE_LIST), POINTER(TREE_LIST), ctypes.c_int, ctypes.c_long]
    lib.best_sos_neighbour.restype = PAIR

    num_trees = len(trees)
    ctreelist = TREE_LIST(num_trees, (TREE * num_trees)(*[t.ctree for t in trees]))
    num_neighbours = len(neighbourhood)
    cneighbourhoud = TREE_LIST(num_neighbours, (TREE * num_neighbours)(*[n.ctree for n in neighbourhood]))

    ret_pair = lib.best_sos_neighbour(cneighbourhoud, ctreelist, n_cores, t_value)

    if ret_pair.index == -1:
        # No better neighbour was found
        raise NoBetterNeighbourFound
    else:
        # Only one neighbour with smallest value found
        return neighbourhood[ret_pair.index], ret_pair.sos

    # print("Multiple choices!")
    # multiple neighbours with same SOS value found, choosing by the given selection method
    # if select == 'first':
    #     return better_neighbours[0], cur_best
    # elif select == 'last':
    #     return better_neighbours[-1], cur_best
    # elif select == 'random':
    #     return choice(better_neighbours), cur_best


def search_neighbourhood_separate(t: TimeTree, trees: TimeTreeSet, t_value: int, n_cores=None, select='first',
                                  rank=False):
    if select not in SELECT_LIST:
        raise ValueError(f"The 'select' parameter should be"
                         f" in {SELECT_LIST} but {select} was given.")

    trees.get_common_clades()

    neighbourhood = intelligent_neighbourhood_separate(t, trees.common_clades, rank)

    better_neighbours = []
    cur_best = t_value
    for n in neighbourhood:
        sos = compute_sos_mt(t=n, trees=trees, n_cores=n_cores)
        if sos < cur_best:
            cur_best = sos
            better_neighbours = [n]
        elif sos == cur_best and sos < t_value:
            better_neighbours.append(n)

    if len(better_neighbours) == 0:
        # No better neighbour was found
        raise NoBetterNeighbourFound
    elif len(better_neighbours) == 1:
        # Only one neighbour with smallest value found
        return better_neighbours[0], cur_best

    # multiple neighbours with same SOS value found, choosing by the given selection method
    if select == 'first':
        return better_neighbours[0], cur_best
    elif select == 'last':
        return better_neighbours[-1], cur_best
    elif select == 'random':
        return choice(better_neighbours), cur_best


def search_hop_neighbourhood_greedy(t, trees, t_value: int, hop_taxa_orders, n_cores=None, select='first'):
    if select not in SELECT_LIST:
        raise ValueError(f"The 'select' parameter should be"
                         f" in {SELECT_LIST} but {select} was given.")

    from TreeVec import TreeVec

    # todo this is just one of the labels, do we need to take all?
    neighbourhood = t[0].hop_neighbourhood()

    better_neighbours = []
    cur_best = t_value

    for n in neighbourhood:
        cur_neighbour_list = []
        # todo need to get the different labels for each neighbour...
        for cur_order in hop_taxa_orders:
            # TreeVec(tree=self.treevec2tree(), leaf2idx=leaf2idx)
            cur_tree = TreeVec(
                tree=n.treevec2tree(),
                leaf2idx=cur_order
            )
            cur_neighbour_list.append(cur_tree)

        sos = compute_hop_sos_mt(t=cur_neighbour_list, trees=trees, n_cores=n_cores)
        if sos < cur_best:
            cur_best = sos
            better_neighbours = [n]
        elif sos == cur_best and sos < t_value:
            better_neighbours.append(n)

    if len(better_neighbours) == 0:
        # No better neighbour was found
        raise NoBetterNeighbourFound
    elif len(better_neighbours) >= 1:
        if len(better_neighbours) > 1:
            print("There was a tie break in neighbours to pick,"
                  " currently not implemented for HOP, picking first in list.")
        better_neighbour_output = []
        for cur_order in hop_taxa_orders:
            cur_tree = TreeVec(
                # newick_str=better_neighbours[0].treevec2newick(),
                tree=better_neighbours[0].treevec2tree(),
                leaf2idx=cur_order
            )
            better_neighbour_output.append(cur_tree)
        return better_neighbour_output, cur_best

    # todo currently not implemented
    # print("Multiple choices!")
    # multiple neighbours with same SOS value found, choosing by the given selection method
    # if select == 'first':
    #     return better_neighbours[0], cur_best
    # elif select == 'last':
    #     return better_neighbours[-1], cur_best
    # elif select == 'random':
    #     return choice(better_neighbours), cur_best