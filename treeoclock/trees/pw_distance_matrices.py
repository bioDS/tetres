__author__ = 'Lena Collienne'
# Computing pairwise distance matrices -- RNNI, RF

import matplotlib.pyplot as plt

from tree_io import *


def pw_rnni_dist(filename):
    tree_list = read_nexus(filename)[0]
    num_trees = tree_list.num_trees

    # Create empty distance matrix as input for pw_distances
    distances = np.zeros(shape=(num_trees,num_trees),dtype=np.int32)
    pw_distances(tree_list, distances)

    # Write computed distance matrix to file -- this can be done way more efficient (at least half the matrix is 0)
    f = open('output/pw_rnni_distance_matrix.txt', 'w+')
    np.savetxt(f, distances, fmt="%i")
    f.close()
    return distances


def pw_rf_dist(filename):
    tree_list,name_dict = read_nexus(filename, ete3=True)
    num_trees = len(tree_list)
    distances = np.zeros(shape=(num_trees,num_trees),dtype=np.int32)
    for i in range(0,num_trees):
        for j in range(i,num_trees):
            distances[i][j] = (tree_list[i].robinson_foulds(tree_list[j])[0])
    f = open('output/pw_rf_distance_matrix.txt', 'w+')
    np.savetxt(f, distances, fmt="%i")
    f.close()
    return distances


def plot_pw_dist_hist(tree_file, dist):
    # Get number of leaves:
    for line in list(open(tree_file)):
        ntax_str = re.search(r'ntax\D*(\d*)', line, re.I)
        if(ntax_str != None):
            num_leaves = int(ntax_str.group(1))
    diameter = int((num_leaves-1)*(num_leaves-2)/2)
    # Get number of trees:
    num_trees = 100
    if (dist == 'rnni'):
        # Compute pairwise distance matrices
        pw_rnni = pw_rnni_dist(tree_file)
        plt.hist(pw_rnni[ pw_rnni!=0 ], bins = int(diameter))
        plt.show()

    elif  (dist == 'rf'):
        pw_rf = pw_rf_dist(tree_file)
        plt.hist(pw_rf[pw_rf!=0], bins = int(diameter))
        plt.show()
