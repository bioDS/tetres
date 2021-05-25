
from ..trees.time_trees import TimeTree, TimeTreeSet
from .compute_sos import compute_sos

def search_neighbourhood_greedy(t: TimeTree, trees: TimeTreeSet, t_value: int, n_cores=None):

    neighbourhood = t.neighbours()
    better_neighbours = []
    for n in neighbourhood:
        sos = compute_sos(t=n, trees=trees, n_cores=n_cores)
        if sos < t_value:
            better_neighbours.append(n)
    # TODO put the neighbours in a priority queue with the sos value ?

    return better_neighbours[1]

