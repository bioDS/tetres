from treeoclock.trees.time_trees import TimeTreeSet, TimeTree
from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary._tree_proposals import search_neighbourhood_greedy, NoBetterNeighbourFound


def greedy(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree):
    # parameter for the starting tree, or just a starting tree given

    initial_sos = compute_sos_mt(start, trees, n_cores=n_cores)
    centroid = start

    while True:
        try:
            centroid, initial_sos = search_neighbourhood_greedy(t=centroid, trees=trees, t_value=initial_sos,
                                                                n_cores=n_cores, select=select)
        except NoBetterNeighbourFound:
            # This is thrown when no neighbour has a better SoS value, i.e. the loop can be stopped
            break

    return centroid
