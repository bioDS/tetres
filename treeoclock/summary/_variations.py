from treeoclock.trees.time_trees import TimeTreeSet, TimeTree
from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary._tree_proposals import search_neighbourhood_greedy, NoBetterNeighbourFound


def greedy(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    sos = compute_sos_mt(start, trees, n_cores=n_cores)
    centroid = start

    while True:
        try:
            centroid, sos = search_neighbourhood_greedy(t=centroid, trees=trees, t_value=sos,
                                                        n_cores=n_cores, select=select)
        except NoBetterNeighbourFound:
            # This is thrown when no neighbour has a better SoS value, i.e. the loop can be stopped
            break
    return centroid, sos


def inc_sub(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):

    sub_size = kwargs["subsample_size"]  # TODO sample size parameter?

    sample = trees  # TODO subsample the set of trees
    sos = compute_sos_mt(start, trees, n_cores=n_cores)

    # subsample.extend(
    #     [copy_trees.pop(rand.randrange(len(copy_trees))) for _ in range(min(subsample_size, len(copy_trees)))])
    # while cp_trees:
    while True:
        # TODO increase subsample
        cen, _ = greedy(trees=sample, n_cores=n_cores, select=select, start=start)
        new_sos = compute_sos_mt(cen, trees, n_cores=n_cores)
        if new_sos <= sos:
            sos = new_sos
        else:
            break


def iter_sub(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    
    sub_size = 200  # TODO sample size parameter?
    
    new_sos = 0
    
    while new_sos <= sos:
        sos = new_sos

        cur_sample = [] # subsample = rand.sample(trees, min(subsample_size, len(trees)))

        new_cen, _ = greedy(cur_sample, n_cores=n_cores, select=select, start=start)
        new_sos = compute_sos_mt(new_cen, trees, n_cores)
    
    
    
    
    return 0
