from multiprocessing import Pool

from treeoclock.trees.time_trees import TimeTree, TimeTreeSet, findpath_distance


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
    d_name = 'RSV2'

    myts = TimeTreeSet(f'/Users/larsberling/Desktop/CodingMA/Git/Summary/MDS_Plots/{d_name}/{d_name}.trees')
    compute_sos(myts[0], myts)
    compute_sos_mp(myts[0], myts)
