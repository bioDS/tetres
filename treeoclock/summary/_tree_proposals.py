
from ..trees.time_trees import TimeTree, TimeTreeSet
from .compute_sos import compute_sos

def search_neighbourhood_greedy(t: TimeTree, trees: TimeTreeSet, t_value: int):

    neighbourhood = t.neighbours()
    better_neighbours = []
    cur_best = t_value
    for n in neighbourhood:
        sos = compute_sos(t=n, trees=trees)
        if sos < cur_best:
            cur_best = sos
            better_neighbours = [n]
        elif sos == cur_best and sos < t_value:
            better_neighbours.append(n)

    if len(better_neighbours) == 1:
        # Only one neighbour with smallest value found
        return better_neighbours[0]

    # TODO select one neighbour with a specified selection = [last, first, random?]
    return better_neighbours
