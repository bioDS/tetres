from multiprocessing import Pool, Array

from treeoclock.trees.time_trees import TimeTree, TimeTreeSet, findpath_distance

global gtrees
gtrees = []


def compute_sos_mp(t: TimeTree, trees: TimeTreeSet, n_cores: int = None) -> int:
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
        dists = p.starmap(findpath_distance, [(t.etree, i.etree) for i in trees])
    return sum(i * i for i in dists)


def compute_sos(t: TimeTree, trees: TimeTreeSet):
    sos = 0
    for i in trees:
        sos += findpath_distance(t.ctree, i.ctree)**2
    return sos

if __name__ == '__main__':
    d_name = 'dispg2d_newzealand_large_time_1'

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')

    from timeit import default_timer as timer

    start = timer()
    f = compute_sos_mp(myts[0], myts)
    print(timer()-start)

    start = timer()
    s = compute_sos(myts[0], myts)
    print(timer()-start)

    print(f, s)
