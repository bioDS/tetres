from tetres.summary.frechet_mean import frechet_mean


def find_convex_hull_trees(chain, i):
    # todo chain index as parameter
    fm = frechet_mean(chain[i].trees)
    dm = chain.pwd_matrix(i)

    fm_d = [fm.fp_distance(t) for t in chain[i].trees]
    convex_list = [True for _ in range(len(chain[i].trees))]
    for t in range(len(chain[i].trees)-1):
        # fix t to be convex, look at all trees closer to fm
        if convex_list[t]:
            for s in range(t,len(chain[i].trees)):
                if convex_list[s]:
                    # if the triangle of t, s and fm has proper lengths than t stays convex, otherwise it is no longer convex
                    if fm_d[t] > fm_d[s] and dm[t, s] < (fm_d[t]/2):
                        convex_list[s] = False
                        # if fm_d[t] > dm[t, s] + fm_d[s]:
    return convex_list
