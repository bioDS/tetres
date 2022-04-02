from treeoclock.trees.time_trees import TimeTreeSet, TimeTree
from treeoclock.summary.compute_sos import compute_sos_mt
from treeoclock.summary._tree_proposals import search_neighbourhood_greedy, NoBetterNeighbourFound, \
    search_neighbourhood_separate, search_neighbourhood_greedy_omp

import random


def greedy(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, tree_log_file='', **kwargs):
    # Chooses always the first tree with better value in the neighbourhood
    sos = compute_sos_mt(start, trees, n_cores=n_cores)
    centroid = start
    count = 0  # Used to name the trees for the logfile
    if tree_log_file:
        with open(tree_log_file, "a") as log:
            log.write(f"tree {count} = {centroid.get_newick()}\n")
            count += 1
    while True:
        try:
            centroid, sos = search_neighbourhood_greedy(t=centroid, trees=trees, t_value=sos,
                                                        n_cores=n_cores, select=select)
        except NoBetterNeighbourFound:
            # This is thrown when no neighbour has a better SoS value, i.e. the loop can be stopped
            break
        if tree_log_file:
            with open(tree_log_file, "a") as log:
                log.write(f"tree {count} = {centroid.get_newick()}\n")
                count += 1

    if tree_log_file:
        with open(tree_log_file, "a") as log:
            log.write("End;")

    return centroid, sos


def greedy_omp(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, tree_log_file='', **kwargs):
    # Chooses always the first tree with better value in the neighbourhood
    sos = compute_sos_mt(start, trees, n_cores=n_cores)
    centroid = start
    count = 0  # Used to name the trees for the logfile
    if tree_log_file:
        with open(tree_log_file, "a") as log:
            log.write(f"tree {count} = {centroid.get_newick()}\n")
            count += 1
    while True:
        try:
            centroid, sos = search_neighbourhood_greedy_omp(t=centroid, trees=trees, t_value=sos,
                                                        n_cores=n_cores, select=select)
        except NoBetterNeighbourFound:
            # This is thrown when no neighbour has a better SoS value, i.e. the loop can be stopped
            break
        if tree_log_file:
            with open(tree_log_file, "a") as log:
                log.write(f"tree {count} = {centroid.get_newick()}\n")
                count += 1

    if tree_log_file:
        with open(tree_log_file, "a") as log:
            log.write("End;")

    return centroid, sos


def inc_sub(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    # Starts with a subset of the trees and then computes the centroid while incrementally increasing the subset

    # Initializing all parameters
    sample = TimeTreeSet()
    centroid = start
    count = 0  # Used to name the trees for the logfile
    if kwargs["tree_log_file"]:
        with open(kwargs["tree_log_file"], "a") as log:
            log.write(f"tree {count} = {centroid.get_newick()}\n")
            count += 1
    sos = compute_sos_mt(start, trees, n_cores=n_cores)

    index_list = list(range(0, len(trees)))
    random.shuffle(index_list)

    iterations = kwargs["max_iterations"]

    while index_list and (iterations is None or iterations > 0):
        if iterations is not None:
            iterations -= 1

        sample.trees.extend(trees.trees[index_list.pop()] for _ in range(min(kwargs["subsample_size"], len(index_list))))

        new_cen, _ = greedy_omp(trees=sample, n_cores=n_cores, select=select, start=centroid)
        new_sos = compute_sos_mt(new_cen, trees, n_cores=n_cores)

        if new_sos <= sos:
            sos = new_sos
            centroid = new_cen
            if kwargs["tree_log_file"]:
                with open(kwargs["tree_log_file"], "a") as log:
                    log.write(f"tree {count} = {centroid.get_newick()}\n")
                    count += 1
        else:
            break

    if kwargs["tree_log_file"]:
        with open(kwargs["tree_log_file"], "a") as log:
            log.write("End;")

    return centroid, sos


def iter_sub(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    # Computes the centroid for a subset and repeats this sampling a new subset but starting computation from
    # previously computed centroid

    # Initializing all parameters
    sample = TimeTreeSet()
    cen = start
    sos = compute_sos_mt(start, trees, n_cores=n_cores)

    index_list = list(range(0, len(trees)))
    random.shuffle(index_list)

    iterations = kwargs["max_iterations"]

    while index_list and (iterations is None or iterations > 0):
        if iterations is not None:
            iterations -= 1

        sample.trees.extend(
            trees.trees[index_list.pop()] for _ in range(min(kwargs["subsample_size"], len(index_list))))

        new_cen, _ = greedy_omp(trees=sample, n_cores=n_cores, select=select, start=cen)
        new_sos = compute_sos_mt(new_cen, trees, n_cores=n_cores)

        if new_sos <= sos:
            sos = new_sos
            cen = new_cen
        else:
            break
    return cen, sos


def separate(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    # Separates the moves and only computes the rank neighbours if the tree contains all clades of the input set

    sos = compute_sos_mt(start, trees, n_cores=n_cores)
    centroid = start

    rank = True
    switch = True
    while True:
        try:
            centroid, sos = search_neighbourhood_separate(t=centroid, trees=trees, t_value=sos,
                                                          n_cores=n_cores, select=select, rank=rank)
            switch = True
        except NoBetterNeighbourFound:
            # This is thrown when no neighbour has a better SoS value, i.e. the loop can be stopped
            if not switch:
                break
            rank = not rank
            switch = False
    return centroid, sos


def onlyone(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    # Prefers to do one move (NNI vs. Rank) and only switches if a local optimum is reached

    sos = compute_sos_mt(start, trees, n_cores=n_cores)
    centroid = start

    # rank = False
    main_move = True  # If False: prefers the NNI move, if True: prefers the Rank move
    switch = True
    while True:
        try:
            # Doing NNI moves
            centroid, sos = search_neighbourhood_separate(t=centroid, trees=trees, t_value=sos,
                                                          n_cores=n_cores, select=select, rank=main_move)
            switch = True
        except NoBetterNeighbourFound:
            # This is thrown when no neighbour has a better SoS value, i.e. the loop can be stopped
            try:
                # Doing one rank move
                centroid, sos = search_neighbourhood_separate(t=centroid, trees=trees, t_value=sos, n_cores=n_cores,
                                                              rank=not main_move)
                switch = True
            except NoBetterNeighbourFound:
                if not switch:
                    break
                # rank = not rank
                switch = False
    return centroid, sos


def update_with_one(trees: TimeTreeSet, n_cores: int, select: str, start: TimeTree, **kwargs):
    # Initializing all parameters
    sample = TimeTreeSet()
    centroid = start

    sos = compute_sos_mt(start, trees, n_cores=n_cores)

    index_list = list(range(0, len(trees)))
    random.shuffle(index_list)

    iterations = kwargs["max_iterations"]
    if iterations is not None:
        iterations -= 1

    sample.trees.extend(trees.trees[index_list.pop()] for _ in range(min(kwargs["subsample_size"], len(index_list))))
    new_cen, _ = greedy_omp(trees=sample, n_cores=n_cores, select=select, start=start)
    new_sos = compute_sos_mt(new_cen, trees, n_cores=n_cores)
    if new_sos < sos:
        sos = new_sos
        centroid = new_cen

    while index_list and (iterations is None or iterations > 0):
        if iterations is not None:
            iterations -= 1

        sample.trees.append(trees.trees[index_list.pop()])  # Adding one tree to the current sample

        new_cen, _ = greedy(trees=sample, n_cores=n_cores, select=select, start=centroid)
        new_sos = compute_sos_mt(new_cen, trees, n_cores=n_cores)

        # d_testing = centroid.fp_distance(new_cen)
        # if d_testing != 0:
        #     print(d_testing)

        if new_sos <= sos:
            sos = new_sos
            centroid = new_cen
        else:
            break

    return centroid, sos
