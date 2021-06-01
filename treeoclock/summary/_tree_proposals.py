from treeoclock.trees.time_trees import TimeTree, TimeTreeSet
from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary._constants import SELECT_LIST

from random import choice


class NoBetterNeighbourFound(Exception):
    """Raised when no neighbour has a better SOS value"""
    pass


def search_neighbourhood_greedy(t: TimeTree, trees: TimeTreeSet, t_value: int, n_cores=None, select='first'):
    neighbourhood = t.neighbours()
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

    raise ValueError(f"The 'select' parameter should be"
                     f" in {SELECT_LIST} but {select} was given.")
