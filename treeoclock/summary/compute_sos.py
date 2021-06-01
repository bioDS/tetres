from multiprocessing.pool import ThreadPool as Pool

from treeoclock.trees.time_trees import TimeTree, TimeTreeSet, findpath_distance


def compute_sos_mt(t: TimeTree, trees: TimeTreeSet, n_cores: int = None) -> int:
    """
    Computes the sum of squared distances for the tree t and the set trees, using n_cores processing cores

    Slower than not using mp as this version needs to compute/convert from etree to ctree again!!!

    :param t: A tree
    :type t: TimeTree
    :param trees: A set of trees
    :type trees: TimeTreeSet
    :param n_cores: Number of processing cores
    :type n_cores: int
    :return: The sum of squared distances
    :rtype: int
    """
    with Pool(n_cores) as p:
        dists = p.starmap(findpath_distance, [(t.ctree, i.ctree) for i in trees])
    return sum(i * i for i in dists)


def compute_sos(t: TimeTree, trees: TimeTreeSet):
    """
        Computes the sum of squared distances for the tree t and the set trees

        Currently faster than using multiple cores, as this uses the ctree structure directly without the need of converting

        :param t: A tree
        :type t: TimeTree
        :param trees: A set of trees
        :type trees: TimeTreeSet
        :return: The sum of squared distances
        :rtype: int
        """
    sos = 0
    for i in trees:
        sos += findpath_distance(t.ctree, i.ctree)**2
    return sos


if __name__ == '__main__':
    d_name = "RSV2"
    
    from timeit import default_timer as timer
    from random import randint

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')
    
    times = []
    times_mp = []
    
    for _ in range(100):
        index = randint(0, len(myts)-1)
        s = timer()
        sos1 = compute_sos(myts[index], myts)
        times.append(timer() - s)

        s = timer()
        sos2 = compute_sos_mt(myts[index], myts)
        times_mp.append(timer() - s)
        
        if not sos1 == sos2:
            print("FAIL")
    
    import numpy as np
    print(np.mean(times), np.mean(times_mp))
    