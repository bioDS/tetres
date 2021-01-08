__author__ = 'Lena Collienne'
# Computing pairwise distance matrices -- RNNI, RF

import ctypes
from tree_io import *
from ete3 import Tree

def pw_rnni_dist(filename):
    # tree_list = read_nexus('../../../rnni_simulations/mcmc_simulation.trees')[0]
    tree_list = read_nexus(filename)[0]
    num_trees = tree_list.num_trees

    # Create empty distance matrix as input for pw_distances
    distance = np.zeros(shape=(num_trees,num_trees),dtype=np.int32)
    pw_distances(tree_list, distance)
    # TODO: convert matrix into numpy matrix / printable format!
    for i in range(0,num_trees-1):
        for j in range(i,num_trees-1):
            print(distance[i][j])
        # print(distance[i][i+1]) # This is probably already a good criterion for convergence? -- can we compare this with RF or other distance measures?
    
    # Write computed distance matrix to file -- this can be done way more efficient (at least half the matrix is 0)
    f = open('test/pw_distance_matrix.txt', 'w+')
    np.savetxt(f, distance, fmt="%i")
    f.close()

def pw_rf_dist(filename):
    tree_list,name_dict = read_nexus(filename, ete3=True)
    for i in range(0,len(tree_list)-1):
        for j in range(i,len(tree_list)-1):
            print(tree_list[i].robinson_foulds(tree_list[j])[0])


if __name__ == '__main__':
    pw_rnni_dist('test/test_9.trees')
    pw_rf_dist('test/test_9.trees')