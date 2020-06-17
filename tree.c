/* author: Lena Collienne
basic algorithms for ranked trees*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

/* trees will be arrays of nodes, ordered according to their ranks (first n nodes are leaves, order there doesn't matter)*/
typedef struct Node{
    int parent;
    int children[2];
} Node;

// List of trees (e.g. as output of NNI move (2 trees) or findpath(d trees))
typedef struct Tree_List{
    int num_trees;
    Node* trees; //?
} Tree_List;


/*read tree  from file*/
Tree_List* read_tree(int num_leaves){

    char filename[100]; // length of this could be variable as well
    printf("What is the file containing trees?\n");
    scanf("%s", filename);
    int num_nodes = num_leaves*2 - 1;

    //TODO: allocate space for the tree_list that can save as many trees as there are lines in this file!
    int num_trees = 2;
    Tree_List * tree_list = malloc(num_trees * sizeof(Tree_List)); // is this correct? 
    for (int i = 0; i < num_trees; i++){
        tree_list[i].trees = malloc(num_nodes * sizeof(Node));
        memset(tree_list[i].trees, -1, num_nodes * sizeof(Node));
        tree_list[i].num_trees = num_trees;
    }

    int *highest_ancestor = malloc(num_leaves * sizeof(int)); // highest_ancestor[i]: index of cluster containing leaf i that is highest below the currently considered cluster
    memset(highest_ancestor, 1, num_leaves * sizeof(int));

    FILE *f;
    f = fopen(filename, "r");

    int max_str_length = 3 * num_leaves * num_nodes; //This is a rough guess for the maximum length of a tree as string. Can we get a better upper bound?

    char *buffer = malloc(max_str_length * sizeof(char));
    memset(buffer, '\0', max_str_length * sizeof(char));

    int current_tree = 0;
    //loop through lines (trees) in file
    while(fgets(buffer, max_str_length * sizeof(char), f) != NULL){
        char *cluster_list, *cluster; //malloc()???
        char * tree_str = malloc(strlen(buffer) * sizeof(char));
        strcpy(tree_str, buffer);
        tree_str[strcspn(tree_str, "\n")] = 0; // delete newline at end of each line that has been read
        int rank = num_leaves;
        //Find clusters
        while((cluster_list = strsep(&tree_str, "}")) != NULL){
            cluster_list += 2; // ignore first two characters [{ or ,{
            if(strlen(cluster_list) > 0){ //ignore last bit (just contained ])
                // Find leaves in clusters
                while((cluster = strsep(&cluster_list, ",")) != NULL){
                    int actual_node = atoi(cluster);
                    // update node relations if current leaf appears for first time
                    if(tree_list[current_tree].trees[actual_node - 1].parent == -1){
                        tree_list[current_tree].trees[actual_node - 1].parent = rank;
                        // update the current internal node (rank) to have the current leaf as child
                        if(tree_list[current_tree].trees[rank].children[0] == -1){
                            tree_list[current_tree].trees[rank].children[0] = actual_node - 1;
                        } else{
                            tree_list[current_tree].trees[rank].children[1] = actual_node - 1;
                        }
                    } else {
                        tree_list[current_tree].trees[highest_ancestor[actual_node - 1]].parent = rank;
                        // update cluster relation if actual_node already has parent assigned (current cluster is union of two clusters or cluster and leaf)
                        if(tree_list[current_tree].trees[rank].children[0] == -1 || tree_list[current_tree].trees[rank].children[0] == highest_ancestor[actual_node - 1]){ // first children should not be assigned yet. I if contains same value, overwrite that one
                            tree_list[current_tree].trees[rank].children[0] = highest_ancestor[actual_node - 1];
                        } else if (tree_list[current_tree].trees[rank].children[1] == -1)
                        {
                            tree_list[current_tree].trees[rank].children[1] = highest_ancestor[actual_node - 1];
                        }
                    }
                    // set new highest ancestor of leaf
                    highest_ancestor[actual_node - 1] = rank;
                }
            }
            rank++;
        }
        current_tree++;
    }

    // //check if read_tree reads trees correctly
    // for (int k = 0; k < num_trees; k++){
    //     for(int i = 0; i < 2 * num_leaves - 1; i++){
    //         if (i < num_leaves){
    //             printf("highest ancestor of node %d has rank %d\n", i, highest_ancestor[i] + 1);
    //             printf("leaf %d has parent %d\n", i+1, tree_list[k].trees[i].parent);
    //         } else{
    //             printf("node %d has children %d and %d\n", i, tree_list[k].trees[i].children[0], tree_list[k].trees[i].children[1]);
    //             printf("leaf %d has parent %d\n", i+1, tree_list[k].trees[i].parent);
    //         }
    //     }
    // }

    fclose(f);
    free(buffer);
    free(highest_ancestor);

    return tree_list;
}


// write tree into given file
void write_tree(Node * tree, int num_leaves, char * filename){

    // TODO: Use malloc!
    char tree_str[5 * num_leaves]; //TODO: replace BUFSIZ with a reasonable upper bound on length of string for tree 

    // create matrix cluster*leaves -- 0 if leaf is not in cluster, 1 if it is in cluster
    int clusters[num_leaves - 1][num_leaves]; // save all clusters as list of lists TODO: malloc
    for (int i = 0; i <num_leaves ; i++){
        for (int j = 0; j < num_leaves - 1; j++){
            clusters[j][i] = 0; //initialise all entries to be 0
        }
        int j = i;
        while (tree[j].parent != -1){
            j = tree[j].parent;
            clusters[j - num_leaves][i] = 1;
        }
        clusters[num_leaves][i] = 1;
    }

    // convert matrix into output string tree_str
    sprintf(tree_str, "[{");
    int tree_str_pos = 1; //last position in tree_str that is filled with a character
    for (int i = 0; i < num_leaves - 1; i++){
        for (int j = 0; j < num_leaves; j++){
            if (clusters[i][j] == 1){
                char leaf_str[100];
                sprintf(leaf_str, "%d,", j + 1);
                // printf("leaf_str: %s\n", leaf_str);
                strcat(tree_str, leaf_str);
                tree_str_pos += 2;
            }
        }
        tree_str[tree_str_pos] = '\0'; // delete last komma
        strcat(tree_str, "},{");
        tree_str_pos +=2;
    }
    tree_str[tree_str_pos] = '\0'; // delete ,{ at end of tree_str
    tree_str[tree_str_pos - 1] = '\0';
    strcat(tree_str, "]");
    printf("%s\n", tree_str);

    // write tree as string to file
    FILE *f;
    f = fopen(filename, "a"); //add tree at end of output file
    fprintf(f, "%s\n", tree_str); //This adds a new line at the end of the file -- we need to be careful when we read from such files!
    fclose(f);
}


int * nni_move(Node * tree, int rank, int num_leaves, int which_child){
    // NNI move on edge bounded by ranks rank - n and rank - n + 1 (n = num_leaves), moving which_child up
    int rank_in_list = rank + num_leaves - 1;
    int num_nodes = 2 * num_leaves - 1;
    if(tree[rank_in_list].parent != rank_in_list + 1){
        printf("Can't do an NNI - interval [%d, %d] is not an edge!\n", rank_in_list, rank_in_list + 1);
    } else{
        int child_moved_up;
        for (int i = 0; i < 2; i++){
            if (tree[rank_in_list+1].children[i] != rank_in_list){ //find the child of the node of rank_in_list k+1 that is not the node of rank_in_list k
                //update parent/children relations to get nni neighbour
                tree[tree[rank_in_list+1].children[i]].parent = rank_in_list; //update parents
                tree[tree[rank_in_list].children[which_child]].parent = rank_in_list+1;
                child_moved_up = tree[rank_in_list].children[which_child];
                tree[rank_in_list].children[which_child] = tree[rank_in_list+1].children[0]; //update children
                tree[rank_in_list+1].children[i] = child_moved_up;
            }
        }
    }
    return 0;
}

Node * rank_move(Node * tree, int rank, int num_leaves){
    // THIS IS STILL LINEAR TIME -- because we need to copy a tree to give two trees as output!
    // Make a rank moves and give the resulting tree as a Tree_List of length 1 (to be consistent with the output of nni_move)
    int num_nodes = 2 * num_leaves - 1;
    int rank_in_list = rank + num_leaves - 1;
    Tree_List * tree_list = malloc(num_nodes * sizeof(Node) + sizeof(int));
    tree_list[0].trees = malloc(num_nodes * sizeof(Node));

    if (tree[rank_in_list].parent == rank_in_list + 1){
        printf("No rank move possible. The interval [%d,%d] is an edge!\n", rank, rank + 1);
    } else{
        // deep copy tree into tree_list
        for (int i = 0; i < num_nodes; i++){
            tree_list[0].trees[i] = tree[i];
        }

        // update nodes in tree
        for (int r = 0; r < 2; r++){
             // updates for node at rank rank + 1 - r
            tree_list[0].trees[rank_in_list + 1 - r].parent = tree[rank_in_list + r].parent;
            tree_list[0].trees[rank_in_list + 1 - r].children[0] = tree[rank_in_list + r].children[0];
            tree_list[0].trees[rank_in_list + 1 - r].children[1] = tree[rank_in_list + r].children[1];

            // update children of node that now have parent with rank rank + 1 - r
            tree_list[0].trees[tree_list[0].trees[rank_in_list + 1 - r].children[0]].parent = rank_in_list + 1 - r;
            tree_list[0].trees[tree_list[0].trees[rank_in_list + 1 - r].children[1]].parent = rank_in_list + 1 - r;

            // update parent of node rank + r
            if (tree[tree[rank_in_list + 1 - r].parent].children[0] == rank_in_list + 1 - r){
                tree_list[0].trees[tree[rank_in_list + 1 - r].parent].children[0] = rank_in_list + r;
            } else{
                tree_list[0].trees[tree[rank_in_list + 1 - r].parent].children[1] = rank_in_list + r;
            }
        }
    }
    return tree_list[0].trees;
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


Node * findpath(Node *start_tree, Node *dest_tree, int num_leaves){
    int current_mrca; //rank of the mrca that needs to be moved down
    Node * current_tree = start_tree; //MALLOC
    for (int i = num_leaves; i < 2 * num_leaves - 1; i++){
        current_mrca = mrca(start_tree, dest_tree[i].children[0], dest_tree[i].children[1]);
        // move current_mrca down
        while(current_mrca != i){
            bool did_nni = false;
            for (int child_index = 0; child_index < 2; child_index++){
                if (did_nni == false && current_tree[current_mrca].children[child_index] == current_mrca - 1){
                    // which child of current_tree[current_mrca].children needs to be moved up? -- that's the one that is not in the current cluster
                    nni_move(current_tree, current_mrca - 1, num_leaves, 1 - child_index);
                    // NNI move on edge [current_mrca - 1 , current_mrca] that decreases the mrca
                    did_nni = true;
                    current_mrca--;
                }
            }
            if (did_nni = false){
                current_tree = rank_move(current_tree, current_mrca - 1, num_leaves);
                current_mrca--;
            }

        }
    }
    write_tree(current_tree, num_leaves, "./output/output.rtree");
}


int main(){
    // TODO: instead of asking for number of leaves, find an upper bound (caterpillar trees)
    int num_leaves;
    printf("How many leaves do your trees have?\n");
    scanf("%d", &num_leaves);

    Tree_List * trees = read_tree(num_leaves);
    int num_nodes;
    num_nodes = 2 * num_leaves - 1;

    remove("./output/output.rtree");
    write_tree(trees[0].trees, num_leaves, "./output/output.rtree"); //TODO: mkdir ./output/ and create file!
    write_tree(trees[1].trees, num_leaves, "./output/output.rtree"); //TODO: mkdir ./output/ and create file!

    // Tree_List * nni_trees = nni_move(trees[1].trees, 1, num_leaves);
    nni_move(trees[1].trees, 2, 5, 0);
    remove("./output/output.rtree");
    write_tree(trees[1].trees, num_leaves, "./output/output.rtree"); //TODO: mkdir ./output/ and create file!
    nni_move(trees[1].trees, 2, 5, 1);
    write_tree(trees[1].trees, num_leaves, "./output/output.rtree"); //TODO: mkdir ./output/ and create file!

    Tree_List * rank_tree = rank_move(trees[0].trees, 1, num_leaves);
    remove("./output/output.rtree");
    write_tree(rank_tree[0].trees, num_leaves, "./output/output.rtree"); //TODO: mkdir ./output/ and create file!
    return 0;
}
