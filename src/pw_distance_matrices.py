__author__ = 'Lena Collienne'
# Computing pairwise distance matrices -- RNNI, RF

import ctypes
from tree_io import *
from ete3 import Tree

def pw_rnni_dist(filename):
    tree_list = read_nexus(filename)[0]
    num_trees = tree_list.num_trees

    # Create empty distance matrix as input for pw_distances
    distance = np.zeros(shape=(num_trees,num_trees),dtype=np.int32)
    pw_distances(tree_list, distance)

    # Write computed distance matrix to file -- this can be done way more efficient (at least half the matrix is 0)
    f = open('output/pw_rnni_distance_matrix.txt', 'w+')
    np.savetxt(f, distance, fmt="%i")
    f.close()


def pw_rf_dist(filename):
    tree_list,name_dict = read_nexus(filename, ete3=True)
    num_trees = len(tree_list)
    print(num_trees)
    distances = np.zeros(shape=(num_trees,num_trees),dtype=np.int32)
    for i in range(0,num_trees):
        for j in range(i,num_trees):
            print(i,j)
            distances[i][j] = (tree_list[i].robinson_foulds(tree_list[j])[0])
    f = open('output/pw_rf_distance_matrix.txt', 'w+')
    np.savetxt(f, distances, fmt="%i")
    f.close()


if __name__ == '__main__':
    tree_file = input("What is the file with trees?\n")
    pw_rnni_dist(tree_file)
    pw_rf_dist(tree_file)