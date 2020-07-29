/* author: Lena Collienne
basic algorithms for ranked trees*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

/* trees will be arrays of nodes, ordered according to their ranks (first n nodes are leaves, order there doesn't matter)*/
typedef struct Node{
    int parent;
    int children[2];
} Node;

// List of trees (e.g. as output of NNI move (2 trees) or findpath(d trees))
typedef struct Tree_List{
    int num_leaves;
    int num_trees;
    Node** trees;
} Tree_List;


int get_num_digits(int integer){
    int n = integer;
    int num_digits = 0;
    while(n != 0){
        n /= 10;
        num_digits++;
    }
    return num_digits;
}

/*read tree  from file*/
Tree_List read_trees(char* filename){

    FILE *f;
    if ((f = fopen(filename, "r"))){
        int num_leaves;
        int num_trees;
        char * buffer0 = malloc(20 * sizeof(char)); // max number of digits for number of leaves and number of trees in file is assumed to be 20
        // read in first two lines: number of leaves and number of trees
        if (fgets(buffer0, 10 * sizeof(char), f) != NULL){
            num_leaves = atoi(buffer0);
            if (num_leaves == 0){
                printf("Error. Determine number of leaves in input file.\n");
            }
        } // TODO: add else to give error when input file does not have the right format
        if (fgets(buffer0, 10 * sizeof(char), f) != NULL){
            num_trees = atoi(buffer0);
            if (num_trees == 0){
                printf("Error. Determine number of trees in input file.\n");
            }
        }
        free(buffer0);
        int num_nodes = num_leaves*2 - 1;
        int num_digits_n = get_num_digits(num_leaves); // number of digits of the int num_leaves
        int max_str_length = 2 * num_leaves * num_leaves * num_digits_n; //upper bound for the maximum length of a tree as string
    
        Tree_List tree_list;
        tree_list.num_leaves = num_leaves;
        tree_list.num_trees = num_trees;
        tree_list.trees = malloc(num_trees*sizeof(Node*));
        for (int i = 0; i < num_trees; i++){
            tree_list.trees[i] = malloc(num_nodes * sizeof(Node));
            for (int j = 0; j < num_nodes; j++){
                tree_list.trees[i][j].parent = -1;
                tree_list.trees[i][j].children[0] = -1;
                tree_list.trees[i][j].children[1] = -1;
            }
            // memset(tree_list.trees[i], -1, num_nodes * sizeof(Node));
        }

        int *highest_ancestor = malloc(num_leaves * sizeof(int)); // highest_ancestor[i]: index of cluster containing leaf i that is highest below the currently considered cluster
        for(int i = 0; i < num_leaves; i++){
            highest_ancestor[i] = 1;
        }
        char *buffer = malloc(max_str_length * sizeof(char*));
        for(int i = 0; i < max_str_length; i++){
            buffer[i] = '\0';
        }
        // memset(buffer, '\0', max_str_length * sizeof(char));
        int current_tree = 0;
        //loop through lines (trees) in file
        while(fgets(buffer, max_str_length * sizeof(char), f) != NULL){
            // Remove white spaces from string buffer
            int l = 0;
            for(int k=0;buffer[k]!='\0';++k)
            {
                if(buffer[k]!=' ')
                    buffer[l++]=buffer[k];
            }
            buffer[l]='\0';

            // allocate memory for strings saving clusters
            char *cluster_list = malloc(max_str_length / (num_leaves - 1)); // maximum length of a cluster as string
            char *tree_str = malloc(max_str_length * sizeof(char));
            char *cluster;
            for(int i = 0; i < max_str_length; i++){
                tree_str[i] = '\0';
            }
            strcpy(tree_str, buffer);
            tree_str[strcspn(tree_str, "\n")] = 0; // delete newline at end of each line that has been read
            int rank = num_leaves;
            //Find clusters
            while((cluster_list = strsep(&tree_str, "}")) != NULL){
                cluster_list += 2; // ignore first two characters [{ or ,{
                if(strlen(cluster_list) > 0){ //ignore last bit (just contained ])
                    // Find leaves in clusters
                    while((cluster = strsep(&cluster_list, ",")) != NULL){
                        int actual_leaf = atoi(cluster);
                        // update node relations if current leaf appears for first time
                        if(tree_list.trees[current_tree][actual_leaf - 1].parent == -1){
                            tree_list.trees[current_tree][actual_leaf - 1].parent = rank;
                            // update the current internal node (rank) to have the current leaf as child
                            if(tree_list.trees[current_tree][rank].children[0] == -1){
                                tree_list.trees[current_tree][rank].children[0] = actual_leaf - 1;
                            } else{
                                tree_list.trees[current_tree][rank].children[1] = actual_leaf - 1;
                            }
                        } else {
                            tree_list.trees[current_tree][highest_ancestor[actual_leaf - 1]].parent = rank;
                            // update cluster relation if actual_leaf already has parent assigned (current cluster is union of two clusters or cluster and leaf)
                            if(tree_list.trees[current_tree][rank].children[0] == -1 || tree_list.trees[current_tree][rank].children[0] == highest_ancestor[actual_leaf - 1]){ // first children should not be assigned yet. I if contains same value, overwrite that one
                                tree_list.trees[current_tree][rank].children[0] = highest_ancestor[actual_leaf - 1];
                            } else if (tree_list.trees[current_tree][rank].children[1] == -1)
                            {
                                tree_list.trees[current_tree][rank].children[1] = highest_ancestor[actual_leaf - 1];
                            }
                        }
                        // set new highest ancestor of leaf
                        highest_ancestor[actual_leaf - 1] = rank;
                    }
                }
                rank++;
            }
            current_tree++;
            free(cluster_list);
            free(tree_str);
        }
        fclose(f);
        free(buffer);
        free(highest_ancestor);

        // //check if read_trees reads trees correctly
        // for (int k = 0; k < num_trees; k++){
        //     for(int i = 9980; i < 10030; i++){
        //         if (i < num_leaves){
        //             // printf("highest ancestor of node %d has rank %d\n", i, highest_ancestor[i] + 1);
        //             printf("leaf %d has parent %d\n", i+1, tree_list.trees[k][i].parent);
        //         } else{
        //             printf("node %d has children %d and %d\n", i, tree_list.trees[k][i].children[0], tree_list.trees[k][i].children[1]);
        //             printf("leaf %d has parent %d\n", i+1, tree_list.trees[k][i].parent);
        //         }
        //     }
        // }

        return tree_list;
    } else{
        printf("Error. File doesn't exist.\n");
    }
}


// write one tree into given file -- runtime quadratic
void write_tree(Node * tree, int num_leaves, char * filename){
    if (tree == NULL){
        printf("Error. Can't write tree. Given tree doesn't exist.\n");
    } else{
        int num_digits_n = get_num_digits(num_leaves); // number of digits of the int num_leaves
        int max_str_length = 2 * num_leaves * num_leaves * num_digits_n; //upper bound for the maximum length of a tree as string
        char *tree_str = malloc(2 * max_str_length * sizeof(char *));

        // create matrix cluster*leaves -- 0 if leaf is not in cluster, 1 if it is in cluster
        int ** clusters = malloc((num_leaves) * sizeof(int *));    
        for (int i = 0; i < num_leaves; i++){
            clusters[i] = malloc((num_leaves - 1) * sizeof(int));
        }

        for (int i = 0; i <num_leaves ; i++){
            for (int j = 0; j < num_leaves - 1; j++){
                clusters[j][i] = 0; //initialise all entries to be 0
            }
            int j = i;
            while (tree[j].parent != -1){
                j = tree[j].parent;
                clusters[j - num_leaves][i] = 1;
            }
            clusters[num_leaves - 1][i] = 1;
        }

        // convert matrix into output string tree_str
        sprintf(tree_str, "[{");
        int tree_str_pos; //last position in tree_str that is filled with a character
        for (int i = 0; i < num_leaves - 1; i++){
            for (int j = 0; j < num_leaves; j++){
                if (clusters[i][j] == 1){
                    char leaf_str[num_digits_n + 1];
                    sprintf(leaf_str, "%d,", j + 1);
                    strcat(tree_str, leaf_str);
                }
            }
            tree_str_pos = strlen(tree_str) - 1;
            tree_str[tree_str_pos] = '\0'; // delete last comma
            strcat(tree_str, "},{");
            tree_str_pos +=2;
        }
        tree_str[tree_str_pos] = '\0'; // delete ,{ at end of tree_str
        tree_str[tree_str_pos - 1] = '\0';
        strcat(tree_str, "]");
        printf("%s\n", tree_str);

        for (int i = 0; i < num_leaves - 1; i++){
            free(clusters[i]);
        }
        // write tree as string to file
        FILE *f;
        f = fopen(filename, "a"); //add tree at end of output file
        fprintf(f, "%s\n", tree_str); //This adds a new line at the end of the file -- we need to be careful when we read from such files!
        fclose(f);
    }
}

void write_trees(Tree_List tree_list, char * filename){
    FILE *f;
    f = fopen(filename, "w");
    fprintf(f, "%d\n", tree_list.num_leaves);
    fprintf(f, "%d\n", tree_list.num_trees);
    fclose(f);
    for (int i = 0; i < tree_list.num_trees; i++){
        write_tree(tree_list.trees[i], tree_list.num_leaves, filename);
    }
}


int nni_move(Node * tree, int rank_in_list, int num_leaves, int child_moves_up){
    // NNI move on edge bounded by rank rank_in_list and rank_in_list + 1, moving child_stays (index) of the lower node up
    if (tree == NULL){
        printf("Error. No RNNI move possible. Given tree doesn't exist.\n");
    } else{
        int num_nodes = 2 * num_leaves - 1;
        if(tree[rank_in_list].parent != rank_in_list + 1){
            printf("Can't do an NNI - interval [%d, %d] is not an edge!\n", rank_in_list, rank_in_list + 1);
            return 1;
        } else{
            int child_moved_up;
            for (int i = 0; i < 2; i++){
                if (tree[rank_in_list+1].children[i] != rank_in_list){ //find the child of the node of rank_in_list k+1 that is not the node of rank_in_list k
                    // update parent/children relations to get nni neighbour
                    tree[tree[rank_in_list+1].children[i]].parent = rank_in_list; //update parents
                    tree[tree[rank_in_list].children[child_moves_up]].parent = rank_in_list+1;
                    child_moved_up = tree[rank_in_list].children[child_moves_up];
                    tree[rank_in_list].children[child_moves_up] = tree[rank_in_list+1].children[i]; //update children
                    tree[rank_in_list+1].children[i] = child_moved_up;
                }
            }
        }
    }
    return 0;
}

int rank_move(Node * tree, int rank_in_list, int num_leaves){
    // Make a rank move on tree between nodes of rank rank and rank + 1 (if possible)
    int num_nodes = 2 * num_leaves - 1;
    if (tree == NULL){
        printf("Error. No rank move possible. Given tree doesn't exist.\n");
        return 1;
    } else{
        if (tree[rank_in_list].parent == rank_in_list + 1){
            printf("Error. No rank move possible. The interval [%d,%d] is an edge!\n", rank_in_list, rank_in_list + 1);
        } else{
            // update parents of nodes that swap ranks
            int upper_parent;
            upper_parent = tree[rank_in_list + 1].parent;
            tree[rank_in_list + 1].parent = tree[rank_in_list].parent;
            tree[rank_in_list].parent = upper_parent;

            for (int i = 0; i < 2; i++){
                // update children of nodes that swap ranks
                int upper_child = tree[rank_in_list + 1].children[i];
                tree[rank_in_list + 1].children[i] = tree[rank_in_list].children[i];
                tree[rank_in_list].children[i] = upper_child;
            }
            for (int i = 0; i < 2; i++){
                // update parents of children of nodes that swap ranks
                tree[tree[rank_in_list + 1].children[i]].parent = rank_in_list + 1; 
                tree[tree[rank_in_list].children[i]].parent = rank_in_list;
            }
            for (int i = 0; i < 2; i ++){
                // update children of parents of nodes that swap rank
                //first case: nodes that swap ranks share a parent. In this case nothing needs to be changed
                if (tree[rank_in_list + 1].parent == tree[rank_in_list].parent){
                    break;
                }
                else{
                    if (tree[tree[rank_in_list + 1].parent].children[i] == rank_in_list){ //parent pointer of tree[rank_in_list + 1] is already set correctly!
                        tree[tree[rank_in_list + 1].parent].children[i] = rank_in_list + 1;
                    }
                    if (tree[tree[rank_in_list].parent].children[i] == rank_in_list + 1){
                        tree[tree[rank_in_list].parent].children[i] = rank_in_list;
                    }
                }
            }
        }
    }
    return 0;
}

int mrca(Node * tree, int node1, int node2){
    // find mrca of nodes with positions node1 and node2 in tree
    int rank1 = node1;
    int rank2 = node2;
    while (rank1 != rank2){
        if (rank1 < rank2){
            rank1 = tree[rank1].parent;
        } else{
            rank2 = tree[rank2].parent;
        }
    }
    return rank1;
}


int ** findpath(Node *start_tree, Node *dest_tree, int num_leaves){
    // returns a path in matrix representation -- explanation in data_structures.md
    int max_dist = ((num_leaves - 1) * (num_leaves - 2))/2 + 1;
    int ** moves = malloc(max_dist * sizeof(int*)); // save moves in a table: each row is move, column 1: rank of lower node bounding the interval of move, column 2: 0,1,2: rank move, nni where children[0] stays, nni where children[1] stays
    for (int i = 0; i < max_dist; i++){
        moves[i] = malloc(2 * sizeof(int));
        moves[i][0] = 0;
        moves[i][1] = 0;
    }
    int path_index = 0; // next position on path that we want to fill with a tree pointer
    if (start_tree == NULL){
        printf("Error. Start tree doesn't exist.\n");
    } else if (dest_tree == NULL){
        printf("Error. Destination tree doesn't exist.\n");
    } else{
        remove("./output/findpath.rtree");
        // write_tree(start_tree, num_leaves, "./output/findpath.rtree"); // this ruins the running time!!!!!!!!
        int current_mrca; //rank of the mrca that needs to be moved down
        Node * current_tree = malloc((2 * num_leaves - 1) * sizeof(Node));
        current_tree =  start_tree;

        for (int i = num_leaves; i < 2 * num_leaves - 1; i++){
            current_mrca = mrca(current_tree, dest_tree[i].children[0], dest_tree[i].children[1]);
            // move current_mrca down
            while(current_mrca != i){
                bool did_nni = false;
                for (int child_index = 0; child_index < 2; child_index++){ // find out if one of the children of current_tree[current_mrca] has rank current_mrca - 1. If this is the case, we want to make an NNI
                    if (did_nni == false && current_tree[current_mrca].children[child_index] == current_mrca - 1){ // do nni if current interval is an edge
                        // check which of the children of current_tree[current_mrca] should move up by the NNI move 
                        bool found_child = false; //indicate if we found the correct child
                        int child_stays; // index of the child of current_tree[current_mrca] that does not move up by an NNI move
                        // find the index of the child of the parent of the node we currently consider -- this will be the index child_stays that we want in the end
                        int current_child_index = dest_tree[i].children[0]; // rank of already existing cluster in both current_tree and dest_tree
                        while (found_child == false){
                            while (current_tree[current_child_index].parent < current_mrca - 1){ // find the x for which dest_tree[i].children[x] is contained in the cluster induced by current_tree[current_mrca - 1]
                                current_child_index = current_tree[current_child_index].parent;
                            }
                            // find the index child_stays
                            if(current_tree[current_child_index].parent == current_mrca - 1){
                                found_child = true;
                                if (current_tree[current_tree[current_child_index].parent].children[0] == current_child_index){
                                    child_stays = 0;
                                } else{
                                    child_stays = 1;
                                }
                            } else{
                                current_child_index = dest_tree[i].children[1];
                            }
                        }
                        nni_move(current_tree, current_mrca - 1, num_leaves, 1 - child_stays);
                        moves[path_index][1] = 1 + child_stays;
                        did_nni = true;
                        current_mrca--;
                    }
                }
                if (did_nni == false){
                    rank_move(current_tree, current_mrca - 1, num_leaves);
                    moves[path_index][1] = 0;
                    current_mrca--;
                }
                moves[path_index][0] = current_mrca;
                path_index++;
            }
        }
    }
    printf("length of path: %d\n", path_index);
    return moves;
}

  
Tree_List return_findpath(Tree_List tree_list){
    // print FP trees into output file, assuming that we want FP between the first two trees in Tree_List trees
    int path_index = 0;
    Node * current_tree = malloc((2 * tree_list.num_leaves - 1) * sizeof(Node)); // deep copy start tree
    for (int i = 0; i < 2 * tree_list.num_leaves - 1; i++){
        current_tree[i] = tree_list.trees[0][i];
    }

    // write_tree(tree_list.trees[1], tree_list.num_leaves, "output/tree.rtree");
    // write_tree(tree_list.trees[0], tree_list.num_leaves, "output/tree.rtree");

    int ** fp = findpath(tree_list.trees[0], tree_list.trees[1], tree_list.num_leaves);
    int diameter = (tree_list.num_leaves - 1) * (tree_list.num_leaves - 2) / 2 + 1; // this is not the diameter, but the number of trees on a path giving the diameter (= diameter + 1)

    // printf("first fp done.\n");

    Tree_List findpath_list; // output: list of trees on FP path
    findpath_list.num_leaves = tree_list.num_leaves;
    findpath_list.trees = malloc(diameter * sizeof(Node *));
    for (int i = 0; i < diameter; i++){
        findpath_list.trees[i] = malloc((2*findpath_list.num_leaves - 1) * sizeof(Node));
    }
    // findpath_list.trees[0] = current_tree;
    for (int i = 0; i < 2 * tree_list.num_leaves - 1; i++){
        findpath_list.trees[0][i] = current_tree[i];
    }

    // create actual path by doing moves starting at tree_list.trees[0] with the information in the matrix returned form fp above
    while(path_index < diameter - 1 && fp[path_index][0] > 0){
        if (fp[path_index][1] == 0){
            rank_move(current_tree, fp[path_index][0], tree_list.num_leaves);
        }
        else if (fp[path_index][1] == 1){
            nni_move(current_tree, fp[path_index][0], tree_list.num_leaves, 1);
        } else{
            nni_move(current_tree, fp[path_index][0], tree_list.num_leaves, 0);
        }
        path_index++;
        // deep copy currently last tree one path
        for (int i = 0; i < 2 * tree_list.num_leaves - 1; i++){
            findpath_list.trees[path_index][i] = current_tree[i];
        }
    }
    findpath_list.num_trees = path_index + 1;
    return findpath_list;
}


int main(){
    char filename[200]; // length of filename set to be 200 char max
    printf("What is the file containing trees?\n");
    scanf("%s", filename);

    Tree_List tree_list = read_trees(filename);
    int num_trees = tree_list.num_trees;
    int num_leaves = tree_list.num_leaves;
    int num_nodes = 2 * num_leaves - 1;

    // // check if read_trees reads trees correctly
    // for (int k = 0; k < num_trees; k++){
    //     for(int i = 0; i < 2 * num_leaves - 1; i++){
    //         if (i < num_leaves){
    //             // printf("highest ancestor of node %d has rank %d\n", i, highest_ancestor[i] + 1);
    //             printf("leaf %d has parent %d\n", i+1, tree_list.trees[k][i].parent);
    //         } else{
    //             printf("node %d has children %d and %d\n", i, tree_list.trees[k][i].children[0], tree_list.trees[k][i].children[1]);
    //             printf("node %d has parent %d\n", i, tree_list.trees[k][i].parent);
    //         }
    //     }
    // }

    // write_trees(tree_list, "./output/output.rtree"); // write given trees into file
    // Tree_List findpath_list = return_findpath(tree_list); // write FP into file
    clock_t start_time = time(NULL);
    int ** fp = findpath(tree_list.trees[0], tree_list.trees[1], tree_list.num_leaves); //run FP
    clock_t end_time = time(NULL);
    printf("Time to compute FP(T,R): %f sec\n", difftime(end_time, start_time));
    // write_trees(findpath_list, "./output/fp.rtree");
    return 0;
}