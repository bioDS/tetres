/*********************************************************************

USAGE GUIDE

For some starting tree T and destination tree R, and some new tree T1 where d(T, T1) == 1.

Precomputation:
	Precomp* precomp_data = precomp(T, R)
O(n^2)
Stores n of the trees along FP(T, R), as well as some information about the movement of nodes (O(n^2) memory).
Also stores d(T, R) in precomp_data->dist.

Finding d(T1, R):
	Update* update_data = find_distance_diff(T, R, move, precomp_data)
Where move is the RNNI edit applied to T to produce T1,
and update_data stores the information required to update the precomputed data structures for following iterations.
The answer is given by update_data->distance_diff = d(T1, R) - d(T, R).

move format (as used by Lena):
	move = [t, x]
	represents an RNNI move on longerval [t, t+1], where:
		x = 0: rank move
		x = 1: NNI move, child at index 0 is moved
		x = 2: NNI move, child at index 1 is moved
See uniform_neighbour function for example.

Following iterations:
To find d(T2, R) where d(T1, T2) == 1, first update the precomputed data structures:
	update_precomp(T, move, precomp_data, update_data)
Then repeat the finding distance step:
	update_data = find_distance_diff(T1, R, move2, precomp_data)
Where move2 is the RNNI edit applied to T1 to produce T2.


**********************************************************************/

/********** CONTENTS **********

GENERAL HELPER FUNCTIONS
void swap(long* a, long* b)
void copy_tree(Tree* T, Tree* R)
void erase(long* arr, long index, long n)
long get_num_digits(long longeger)

RNNI HELPER FUNCTIONS
long nni_move(Tree * input_tree, long rank_in_list, long child_moves_up)
long rank_move(Tree * input_tree, long rank_in_list)
void perform_move(Tree* T, long* move)
long mrca(Tree * input_tree, long node1, long node2)
long* uniform_neighbour(Tree * input_tree)

ONLINE RNNI HELPER FUNCTIONS
long find_head(long node, long edit_rank, Tree* T)
long get_members(long head, long cluster, Tree* T, long* members, long n_members)
void label_leaves(long head, Tree* T, char* labels, char lbl)
void label_members(long head, long cluster, Tree* T, char* labels, char lbl)
long find_label(long node, char* labels, Tree* T)
bool contains_label(long node, Tree* T, char* labels, char lbl)
long have_swapped(long edit_node, long cluster, Precomp* precomp_data)

ONLINE RNNI PRECOMPUTATION FUNCTIONS
Precomp* initialise_precomp(long n)
Precomp* precomp(Tree *start_tree, Tree *dest_tree)
Update* initialise_update(long n)
void update_precomp(Tree* T_next, long* move, Precomp* precomp_data, Update* update_data)

ONLINE RNNI DISTANCE COMPUTATION FUNCTIONS
void get_distance_diff_NNI(Tree* T0, Tree* T1, Tree* R, char* labels, long edit_node, Precomp* precomp_data, Update* update_data)
void get_distance_diff_rank(Tree* T0, Tree* T1, Tree* R, char* labels, long edit_node, Precomp* precomp_data, Update* update_data)
Update* find_distance_diff(Tree* T0, Tree* R, long* move, Precomp* precomp_data)

TESTING FUNCTIONS
void generate_random_tree(long n, Tree* T)
char* tree_to_string(Tree * input_tree)
long findpath_distance(Tree *start_tree, Tree *dest_tree)
void test_random_single_edit(long n, long trials)
void test_random_single_edit_chain(long n, long chain_length, long trials)

********** CONTENTS END **********/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <assert.h>

const long DIST_DIFF = 0, PATHS_JOIN = 1, CHILD_MOVES = 2; // indices

typedef struct Node {
	long parent;
	long children[2];
	long time;
} Node;

typedef struct Tree{
	long num_leaves;
	Node * tree;
	long root_time;
} Tree;

typedef struct Path{
	long ** moves;
	long length;
} Path;

typedef struct Precomp{
	Tree** tree_changes;
	long** rank_changes;
	long** rank_identities;
	long dist;
} Precomp;

typedef struct Update{
	long* swaps;
	long distance_diff;
	long paths_join;
	long child_moves;
} Update;


// ********** GENERAL HELPER FUNCTIONS **********

void swap(long* a, long* b) {
	long temp = *a;
	*a = *b;
	*b = temp;
}

void copy_tree(Tree* T, Tree* R) {
	long n = T->num_leaves;
	R->num_leaves = n;
	R->tree = (Node*)malloc((2 * n - 1) * sizeof(Node));
    for (long i = 0; i < 2 * n - 1; i++){
        R->tree[i] = T->tree[i];
    }
}

void erase(long* arr, long index, long n) {
	for (long i = index; i < n-1; i++) {
		arr[i] = arr[i+1];
	}
}

// Number of digits of an longeger -- needed to get an upper bound of the length of an input tree as string (when reading from a file)
long get_num_digits(long longeger) {
	// Lena original
    long n = longeger;
    long num_digits = 0;
    while(n != 0){
        n /= 10;
        num_digits++;
    }
    return num_digits;
}

// ********** GENERAL HELPER FUNCTIONS END **********

// ********** RNNI HELPER FUNCTIONS **********

long nni_move(Tree * input_tree, long rank_in_list, long child_moves_up) {
	// Lena original
    if (input_tree->tree == NULL){
        printf("Error. No RNNI move possible. Given tree doesn't exist.\n");
    } else{
        if(input_tree->tree[rank_in_list].parent != rank_in_list + 1){
            printf("Can't do an NNI - longerval [%ld, %ld] is not an edge!\n", rank_in_list+1, rank_in_list+1 + 1);
            return 1;
        } else{
            long child_moved_up;
            for (long i = 0; i < 2; i++){
                if (input_tree->tree[rank_in_list+1].children[i] != rank_in_list){ //find the child of the node of rank_in_list k+1 that is not the node of rank_in_list k
                    // update parent/children relations to get nni neighbour
                    input_tree->tree[input_tree->tree[rank_in_list+1].children[i]].parent = rank_in_list; //update parents
                    input_tree->tree[input_tree->tree[rank_in_list].children[child_moves_up]].parent = rank_in_list+1;
                    child_moved_up = input_tree->tree[rank_in_list].children[child_moves_up];
                    input_tree->tree[rank_in_list].children[child_moves_up] = input_tree->tree[rank_in_list+1].children[i]; //update children
                    input_tree->tree[rank_in_list+1].children[i] = child_moved_up;
                    break;
                }
            }
        }
    }
    return 0;
}

// Make a rank move on tree between nodes of rank rank and rank + 1 (if possible)
long rank_move(Tree * input_tree, long rank_in_list){
    // Lena original
    if (input_tree->tree == NULL){
        printf("Error. No rank move possible. Given tree doesn't exist.\n");
        return 1;
    } else{
        if (input_tree->tree[rank_in_list].parent == rank_in_list + 1){
            printf("Error. No rank move possible. The longerval [%ld,%ld] is an edge!\n", rank_in_list, rank_in_list + 1);
        } else{
            // update parents of nodes that swap ranks
            long upper_parent;
            upper_parent = input_tree->tree[rank_in_list + 1].parent;
            input_tree->tree[rank_in_list + 1].parent = input_tree->tree[rank_in_list].parent;
            input_tree->tree[rank_in_list].parent = upper_parent;

            for (long i = 0; i < 2; i++){
                // update children of nodes that swap ranks
                long upper_child = input_tree->tree[rank_in_list + 1].children[i];
                input_tree->tree[rank_in_list + 1].children[i] = input_tree->tree[rank_in_list].children[i];
                input_tree->tree[rank_in_list].children[i] = upper_child;
            }
            for (long i = 0; i < 2; i++){
                // update parents of children of nodes that swap ranks
                input_tree->tree[input_tree->tree[rank_in_list + 1].children[i]].parent = rank_in_list + 1; 
                input_tree->tree[input_tree->tree[rank_in_list].children[i]].parent = rank_in_list;
            }
            for (long i = 0; i < 2; i ++){
                // update children of parents of nodes that swap rank
                //first case: nodes that swap ranks share a parent. In this case nothing needs to be changed
                if (input_tree->tree[rank_in_list + 1].parent == input_tree->tree[rank_in_list].parent){
                    break;
                }
                else{
                    if (input_tree->tree[input_tree->tree[rank_in_list + 1].parent].children[i] == rank_in_list){ //parent polonger of input_tree->tree[rank_in_list + 1] is already set correctly!
                        input_tree->tree[input_tree->tree[rank_in_list + 1].parent].children[i] = rank_in_list + 1;
                    }
                    if (input_tree->tree[input_tree->tree[rank_in_list].parent].children[i] == rank_in_list + 1){
                        input_tree->tree[input_tree->tree[rank_in_list].parent].children[i] = rank_in_list;
                    }
                }
            }
        }
    }
    return 0;
}

void perform_move(Tree* T, long* move) {
	if (move[1] == 0) {
		rank_move(T, move[0]);
	} else if (move[1] == 1) {
		nni_move(T, move[0], 0);
	} else {
		nni_move(T, move[0], 1);
	}
}

// find mrca of nodes with positions node1 and node2 in tree
long mrca(Tree * input_tree, long node1, long node2){
	// Lena original
    long rank1 = node1;
    long rank2 = node2;
    while (rank1 != rank2){
        if (rank1 < rank2){
            rank1 = input_tree->tree[rank1].parent;
        } else{
            rank2 = input_tree->tree[rank2].parent;
        }
    }
    return rank1;
}

long* uniform_neighbour(Tree * input_tree){
	// Lena original
    // Perform a random RNNI move (at uniform) on input_tree
    // Deep copy input tree, so we can perform move on it
    long num_leaves = input_tree->num_leaves;
    // Tree * output_tree = malloc(sizeof(Node*) + 3 * sizeof(long));
    // output_tree->num_leaves = num_leaves;
    // output_tree->tree = malloc((2 * num_leaves - 1) * sizeof(Node)); // deep copy start tree
    // for (long i = 0; i < 2 * num_leaves - 1; i++){
    //     output_tree->tree[i] = input_tree->tree[i];
    // }
    // Count number of possible moves (rank longerval + 2*NNI longerval)
    long num_moves = 0; // total number of possible moves on given tree (2* #edge longervals + 1 * #rank longervals)
    long ** move_list = (long**)malloc(2 * (num_leaves - 1) * sizeof(long*));
    for (long i = 0; i < 2*(num_leaves - 1); i++){ // max number of moves is reached if every longernal edge has length one (caterpillar)
        move_list[i] = (long*)malloc(2*sizeof(long));
        move_list[i][0] = -1; // lower node of edge for move
        move_list[i][1] = -1; // rank vs nni move
    }
    // Fill move list
    for (long i = num_leaves; i < 2 * num_leaves - 1; i++){
        // printf("rank: %ld, children[0]: %ld, children[1]: %ld, parent: %ld\n", i, input_tree->tree[i].children[0], input_tree->tree[i].children[1], input_tree->tree[i].parent);
        if (input_tree->tree[i].parent == i+1){
            // printf("move_list index: %ld, num_moves: %ld, type of move: %ld\n", i, num_moves, 1);
            move_list[num_moves][0] = i;
            move_list[num_moves][1] = 1; // NNI move 0
            move_list[num_moves + 1][0] = i;
            move_list[num_moves + 1][1] = 2; // NNI move 1
            // printf("move_list index: %ld, num_moves: %ld, type of move: %ld\n", i, num_moves+1, 2);
            num_moves += 2;
        } else{
            // printf("move_list index: %ld, num_moves: %ld, type of move: %ld\n", i, num_moves, 0);
            move_list[num_moves][0] = i;
            move_list[num_moves][1] = 0; // rank move is 0
            num_moves += 1;
        }
        // printf("num_moves: %ld\n", num_moves);
    }

    // Pick random move
    long r = rand() % (num_moves-1);
    // printf("r: %ld\n", r);
    if (move_list[r][1] == 0){
        rank_move(input_tree, move_list[r][0]);
    } else if (move_list[r][1] == 1){
        nni_move(input_tree, move_list[r][0], 0);
    } else{
        nni_move(input_tree, move_list[r][0], 1);
    }
    long* move = (long*)malloc(2*sizeof(long));
    move[0] = move_list[r][0];
    move[1] = move_list[r][1];
    // free move_list
    for (long i = 0; i < 2*(num_leaves - 1); i++){ // max number of moves is reached if every longernal edge has length one (caterpillar)
        free(move_list[i]);
    }
    free(move_list);
    return move;
}

// ********** RNNI HELPER FUNCTIONS END **********

// ********** ONLINE RNNI HELPER FUNCTIONS **********

long find_head(long node, long edit_rank, Tree* T) {
	// O(n) worst case, overall O(n) total even if called n times
	assert(node < edit_rank);
	long head = node;
	while (T->tree[head].parent <= edit_rank) {
		head = T->tree[head].parent;
	}
	assert(0 <= head && head < 2*T->num_leaves);
	return head;
}

long get_members(long head, long cluster, Tree* T, long* members, long n_members) {
	// O(n) worst case, overall O(n) total even if called n times
	assert(head != -1);
	if (head <= cluster+T->num_leaves) {
		members[n_members] = head;
		n_members++;
	}
	if (head >= T->num_leaves) {
		n_members = get_members(T->tree[head].children[0], cluster, T, members, n_members);
		n_members = get_members(T->tree[head].children[1], cluster, T, members, n_members);
	}
	return n_members;
}

void label_leaves(long head, Tree* T, char* labels, char lbl) {
	// O(n) worst case, overall O(n) total even if called n times
	if (head < T->num_leaves) {
		labels[head] = lbl;
		return;
	}
	assert(T->tree[head].children[0] != -1);
	label_leaves(T->tree[head].children[0], T, labels, lbl);
	label_leaves(T->tree[head].children[1], T, labels, lbl);
}

void label_members(long head, long cluster, Tree* T, char* labels, char lbl) {
	// O(n) worst case, overall O(n) total even if called n times
	long n = T->num_leaves;
	long* members = (long*)malloc((2 * n - 1) * sizeof(long));
	long n_members = get_members(head, cluster, T, members, 0);
	for (long i = 0; i < n_members; i++) {
		labels[members[i]] = lbl;
	}
	labels[head] = lbl;
}

long find_label(long node, char* labels, Tree* T) {
	// O(1)
	if (node < T->num_leaves) {
		return labels[node];
	}
	assert(labels[T->tree[node].children[0]] == labels[T->tree[node].children[1]]);
	return labels[node] = labels[T->tree[node].children[0]];
}

bool contains_label(long node, Tree* T, char* labels, char lbl) {
	// TODO complexity? usage?
	if (node < T->num_leaves) {
		return labels[node] == lbl;
	}
	return contains_label(T->tree[node].children[0], T, labels, lbl) || contains_label(T->tree[node].children[1], T, labels, lbl);
}

long have_swapped(long edit_node, long cluster, Precomp* precomp_data) {
	// TODO this is confusing and convoluted
	// O(1)
	Tree** tree_changes = precomp_data->tree_changes;
	long** rank_changes = precomp_data->rank_changes;
	long** rank_identities = precomp_data->rank_identities;
	if (cluster == 0) return 0;
	//long n = tree_changes[0]->num_leaves;
	long count_same = 0;
	for (long i = 0; i < 2; i++) {
		long edit_rank = rank_changes[cluster-1][edit_node];
		long prev_id = rank_identities[cluster-1][tree_changes[cluster-1]->tree[edit_rank].children[i]];
		//long prev_id = find_identity(tree_changes[cluster-1]->tree[edit_rank].children[i], cluster-1, rank_changes, n);
		edit_rank = rank_changes[cluster][edit_node];
		long current_id = rank_identities[cluster][tree_changes[cluster]->tree[edit_rank].children[i]];
		//long current_id = find_identity(tree_changes[cluster]->tree[edit_rank].children[i], cluster, rank_changes, n);
		if (prev_id == current_id) count_same++;
	}
	if (count_same > 0) {
		return 0;
	} else {
		return 1;
	}
}

// ********** ONLINE RNNI HELPER FUNCTIONS END **********

// ********** ONLINE RNNI PRECOMPUTATION FUNCTIONS **********

Precomp* initialise_precomp(long n) {
	long** rank_identities = (long**)malloc(n * sizeof(long*));
	long** rank_changes = (long**)malloc(n * sizeof(long*));
	Tree** tree_changes = (Tree**)malloc(n * sizeof(Tree*));
	for (long i = 0; i < n; i++) {
		rank_identities[i] = (long*)malloc((n * 2 - 1) * sizeof(long));
		rank_changes[i] = (long*)malloc((n * 2 - 1) * sizeof(long));
		tree_changes[i] = (Tree*)malloc(sizeof(Tree));
		for (long j = 0; j < n * 2 - 1; j++) {
			rank_identities[i][j] = j;
			rank_changes[i][j] = j;
		}
	}
	Precomp* precomp_data = (Precomp*)malloc(sizeof(Precomp));
	precomp_data->tree_changes = tree_changes;
	precomp_data->rank_changes = rank_changes;
	precomp_data->rank_identities = rank_identities;
	return precomp_data;
}

Precomp* precomp(Tree *start_tree, Tree *dest_tree) {
	// Lena's findpath function with some edits to store required data
	// tree_changes: stores first tree in each iteration of findpath (n of the trees along the path)
	// rank_changes:
	// rank_identities:
    long num_leaves = start_tree->num_leaves;
    Precomp* precomp_data = initialise_precomp(num_leaves);
    Tree** tree_changes = precomp_data->tree_changes;
    long** rank_changes = precomp_data->rank_changes;
    long** rank_identities = precomp_data->rank_identities;
    //long max_dist = ((num_leaves - 1) * (num_leaves - 2))/2 + 1;

    //printf();

    long current_mrca; //rank of the mrca that needs to be moved down
    Tree current_tree;
    current_tree.tree = (Node*)malloc((2 * num_leaves - 1) * sizeof(Node));
    current_tree.num_leaves = num_leaves;
    for (long i = 0; i < 2 * num_leaves - 1; i++){
        current_tree.tree[i] = start_tree->tree[i];
    }
    Tree * current_tree_polonger;
    current_tree_polonger = &current_tree;

    long distance = 0;
    for (long i = num_leaves; i < 2 * num_leaves - 1; i++){
    	// i-th cluster of destination tree
    	// ***** PRECOMP *****
    	// TODO this tree_changes should be freed at some polong?
    	copy_tree(current_tree_polonger, tree_changes[i-num_leaves]);
    	for (long j = 0; j < num_leaves*2-1; j++) {
    		rank_identities[i-num_leaves+1][j] = rank_identities[i-num_leaves][j];
    	}
    	// ***************
    	
        current_mrca = mrca(current_tree_polonger, dest_tree->tree[i].children[0], dest_tree->tree[i].children[1]);

        // move current_mrca down
		while(current_mrca != i){
			distance += 1;
            bool did_nni = false;
            for (long child_index = 0; child_index < 2; child_index++){
                if (did_nni == false && current_tree.tree[current_mrca].children[child_index] == current_mrca - 1){
                    bool found_child = false;
                    long child_stays;
                    long current_child_index = dest_tree->tree[i].children[0];
                    while (found_child == false){
                        while (current_tree.tree[current_child_index].parent < current_mrca - 1){
                            current_child_index = current_tree.tree[current_child_index].parent;
                        }
                        // find the index child_stays
                        if(current_tree.tree[current_child_index].parent == current_mrca - 1){
                            found_child = true;
                            if (current_tree.tree[current_tree.tree[current_child_index].parent].children[0] == current_child_index){
                                child_stays = 0;
                            } else{
                                child_stays = 1;
                            }
                        } else{
                            current_child_index = dest_tree->tree[i].children[1];
                        }
                    }
                    nni_move(current_tree_polonger, current_mrca - 1, 1 - child_stays);
                    did_nni = true;
                    
                    assert(current_mrca > 0);
                	assert(current_mrca < 2*num_leaves-1);
                	// ***** PRECOMP *****
                    swap(&rank_identities[i-num_leaves+1][current_mrca-1],
                    		&rank_identities[i-num_leaves+1][current_mrca]);
                    // ***************
                    
                    current_mrca--;
                }
            }
            if (did_nni == false){
                rank_move(current_tree_polonger, current_mrca - 1);
                
                assert(current_mrca > 0);
                assert(current_mrca < 2*num_leaves-1);
                // ***** PRECOMP *****
                swap(&rank_identities[i-num_leaves+1][current_mrca-1],
                    	&rank_identities[i-num_leaves+1][current_mrca]);
                // ***************
                
                current_mrca--;
            }
		}
    }
    free(current_tree.tree);
    // ***** PRECOMP *****
    for (long i = 0; i < num_leaves-1; i++) {
    	for (long j = 0; j < num_leaves*2-1; j++) {
    		rank_changes[i][rank_identities[i][j]] = j;
    	}
    }
    precomp_data->dist = distance;
    // ***************
    return precomp_data;
}

Update* initialise_update(long n) {
	long* swaps = (long*)malloc((n-1) * sizeof(long));
	for (long i = 0; i < n-1; i++) swaps[i] = -1;
	Update* update_data = (Update*)malloc(sizeof(Update));
	update_data->swaps = swaps;
	update_data->distance_diff = 2;
	update_data->paths_join = -1;
	update_data->child_moves = -1;
	return update_data;
}

void update_precomp(Tree* T_next, long* move, Precomp* precomp_data, Update* update_data) {
	long n = T_next->num_leaves;
	Tree** tree_changes = precomp_data->tree_changes;
	long** rank_changes = precomp_data->rank_changes;
	long** rank_identities = precomp_data->rank_identities;
	long* next_move = (long*)malloc(2 * sizeof(long));
	next_move[1] = move[1];
	long edit_node = move[0];
	for (long i = 0; i < n-1; i++) {
		if (i <= update_data->paths_join) {
			next_move[0] = rank_changes[i][edit_node];
			if (next_move[1] != 0) {
				assert(update_data->swaps[i] != -1);
				if (update_data->swaps[i] == 1) {
					next_move[1] = 3 - next_move[1];
				}
			} else if (tree_changes[i]->tree[next_move[0]].parent == next_move[0] + 1) {
				assert(update_data->child_moves != -1);
				next_move[1] = (1-update_data->child_moves) + 1;
			}
			perform_move(tree_changes[i], next_move);
		} else if (update_data->distance_diff != 0) {
			swap(&rank_identities[i][rank_changes[i][edit_node]], &rank_identities[i][rank_changes[i][edit_node+1]]);
			swap(&rank_changes[i][edit_node], &rank_changes[i][edit_node+1]);
		}
	}
}

// ********** ONLINE RNNI PRECOMPUTATION FUNCTIONS END **********

// ********** ONLINE RNNI DISTANCE COMPUTATION FUNCTIONS **********

void get_distance_diff_NNI(Tree* T0, Tree* T1, Tree* R, char* labels, long edit_node, Precomp* precomp_data, Update* update_data) {
	long n = T0->num_leaves;
	Tree** tree_changes = precomp_data->tree_changes;
	long** rank_changes = precomp_data->rank_changes;
	//long** rank_identities = precomp_data->rank_identities;
	long edit_rank;
	for (long i = n; i < 2*n-1; i++) {
		edit_rank = rank_changes[i-n][edit_node];
		long c0 = R->tree[i].children[0];
		long c1 = R->tree[i].children[1];
		labels[c0] = find_label(c0, labels, R);
		labels[c1] = find_label(c1, labels, R);
		
		update_data->paths_join = i-n;
		update_data->swaps[i-n] = have_swapped(edit_node, i-n, precomp_data);
		
		// CASES
		if (labels[c0] > labels[c1]) {
			swap(&c0, &c1);
		}
		if (labels[c0] == labels[c1]) {
			labels[i] = labels[c0];
		} else if (labels[c0] == 'A' && labels[c1] == 'B') {
			update_data->distance_diff = 1;
			return;
		} else if (labels[c0] == 'B' && labels[c1] == 'C') {
			update_data->distance_diff = -1;
			return;
		} else if (labels[c0] == 'A' && labels[c1] == 'C') {
			update_data->distance_diff = 0;
			return;
		} else {
			long head = find_head(c1, edit_rank, tree_changes[i-n]);
			label_members(head, i-n, tree_changes[i-n], labels, labels[c0]);
			labels[i] = labels[c0];
		}
	}
	printf("something has gone terribly wrong\n");
}

void get_distance_diff_rank(Tree* T0, Tree* T1, Tree* R, char* labels, long edit_node, Precomp* precomp_data, Update* update_data) {
	long n = T0->num_leaves;
	Tree** tree_changes = precomp_data->tree_changes;
	long** rank_changes = precomp_data->rank_changes;
	// long** rank_identities = precomp_data->rank_identities;
	long edit_rank = edit_node;
	long child_moves = -1;
	for (long i = n; i < 2*n-1; i++) {
		edit_rank = rank_changes[i-n][edit_node];
		long c0 = R->tree[i].children[0];
		long c1 = R->tree[i].children[1];
		labels[c0] = find_label(c0, labels, R);
		labels[c1] = find_label(c1, labels, R);
		
		update_data->paths_join = i-n;
		
		// CASES
		if (labels[c0] > labels[c1]) {
			swap(&c0, &c1);
		}
		if (labels[c0] == labels[c1]) {
			labels[i] = labels[c0];
		} else if (labels[c0] == 'A' && labels[c1] == 'B') {
			update_data->distance_diff = 1;
			return;
		} else if (labels[c0] == 'C' && labels[c1] == 'D') {
			update_data->distance_diff = -1;
			return;
		} else if (labels[c1] == 'E') {
			assert(i < 2*n-1 - 1);
			long head = find_head(c1, edit_rank, tree_changes[i-n]);
			label_members(head, i-n, tree_changes[i-n], labels, labels[c0]);
			labels[i] = labels[c0];
		} else {
			assert(i < 2*n-1 - 1);
			Tree* next_T1 = (Tree*)malloc(sizeof(Tree));
			Tree* next_T0 = tree_changes[i-n+1];
			copy_tree(next_T0, next_T1);
			long next_edit_rank = rank_changes[i-n+1][edit_node];
			
			char lbl = 'A';
			if (labels[c0] == 'A') {
				lbl = 'B';
			}
			long moving_child = next_T0->tree[next_edit_rank].children[0];
			bool move_zero = !contains_label(moving_child, next_T0, labels, lbl);
			
			for (long i = 0; i < 2*n-1; i++) {
				labels[i] = 'D';
			}
			long head = next_T0->tree[next_edit_rank].parent;
			if (next_T0->tree[head].children[0] == next_edit_rank) {
				head = next_T0->tree[head].children[1];
			} else {
				head = next_T0->tree[head].children[0];
			}
			label_leaves(head, next_T0, labels, 'C');
			if (move_zero) {
				child_moves = 0;
				nni_move(next_T1, next_edit_rank, 0);
				head = next_T0->tree[next_edit_rank].children[0];
				label_leaves(head, next_T0, labels, 'B');
				head = next_T0->tree[next_edit_rank].children[1];
				label_leaves(head, next_T0, labels, 'A');
			} else {
				child_moves = 1;
				nni_move(next_T1, next_edit_rank, 1);
				head = next_T0->tree[next_edit_rank].children[0];
				label_leaves(head, next_T0, labels, 'A');
				head = next_T0->tree[next_edit_rank].children[1];
				label_leaves(head, next_T0, labels, 'B');
			}
			update_data->child_moves = child_moves;
			get_distance_diff_NNI(next_T0, next_T1, R, labels, edit_node, precomp_data, update_data);
			return;
		}
	}
	printf("something has gone terribly wrong\n");
}

Update* find_distance_diff(Tree* T0, Tree* R, long* move, Precomp* precomp_data) {
	long n = T0->num_leaves;
	Update* update_data = initialise_update(n);
	Tree* T1 = (Tree*)malloc(sizeof(Tree));
	copy_tree(T0, T1);
	perform_move(T1, move);
	char* labels = (char*)malloc((2 * n - 1) * sizeof(char));
	for (long i = 0; i < 2*n-1; i++) {
		labels[i] = 'E';
	}
	if (move[1] == 0) {
		label_leaves(T0->tree[move[0]].children[0], T0, labels, 'A');
		label_leaves(T0->tree[move[0]].children[1], T0, labels, 'B');
		label_leaves(T0->tree[move[0]+1].children[0], T0, labels, 'C');
		label_leaves(T0->tree[move[0]+1].children[1], T0, labels, 'D');
		get_distance_diff_rank(T0, T1, R, labels, move[0], precomp_data, update_data);
	} else {
		if (move[1] == 1) {
			label_leaves(T0->tree[move[0]].children[0], T0, labels, 'A');
			label_leaves(T0->tree[move[0]].children[1], T0, labels, 'B');
		} else {
			label_leaves(T0->tree[move[0]].children[0], T0, labels, 'B');
			label_leaves(T0->tree[move[0]].children[1], T0, labels, 'A');
		}
		long sibling = T0->tree[move[0]].parent;
		if (T0->tree[sibling].children[0] == move[0]) {
			sibling = T0->tree[sibling].children[1];
		} else {
			sibling = T0->tree[sibling].children[0];
		}
		label_leaves(sibling, T0, labels, 'C');
		get_distance_diff_NNI(T0, T1, R, labels, move[0], precomp_data, update_data);
	}
	return update_data;
}

// ********** ONLINE RNNI DISTANCE COMPUTATION FUNCTIONS END **********

// ********** TESTING FUNCTIONS **********
//void generate_random_tree(long n, Tree* T) {
//	T->num_leaves = n;
//	T->tree = (Node*)malloc((2 * n - 1) * sizeof(Node));
//	long* avail = (long*)malloc(n * sizeof(long));
//	long n_avail = 0;
//	for (long i = 0; i < n; i++) {
//		avail[n_avail] = i;
//		n_avail++;
//		T->tree[i].children[0] = -1;
//		T->tree[i].children[1] = -1;
//	}
//	T->tree[2*n-2].parent = -1;
//	for (long i = n; i < 2*n-1; i++) {
//    	long r = rand() % n_avail;
//    	long c0 = avail[r];
//    	erase(avail, r, n_avail);
//    	n_avail--;
//    	r = rand() % n_avail;
//    	long c1 = avail[r];
//    	erase(avail, r, n_avail);
//    	n_avail--;
//    	avail[n_avail] = i;
//    	n_avail++;
//
//    	T->tree[i].children[0] = c0;
//    	T->tree[i].children[1] = c1;
//    	T->tree[c0].parent = i;
//    	T->tree[c1].parent = i;
//	}
//	assert(n_avail == 1);
//}

// Return tree as string in cluster format -- for testing purposes
//char* tree_to_string(Tree * input_tree){
//	// Lena original
//    if (input_tree->tree == NULL){
//        printf("Error. Can't write tree. Given tree doesn't exist.\n");
//        return(NULL);
//    } else{
//        long num_leaves = input_tree->num_leaves;
//        long num_digits_n = get_num_digits(input_tree->num_leaves); // number of digits of the long num_leaves
//        long max_str_length = 2 * num_leaves * num_leaves * num_digits_n; //upper bound for the maximum length of a tree as string
//        char *tree_str = (char*)malloc(2 * max_str_length * sizeof(char));
//
//        // Check if input tree is 'correct'
//        // for (long i = 0; i < 2 * num_leaves - 1; i++){
//        //     printf("Node %ld, Parent %ld, Children %ld and %ld\n", i, input_tree->tree[i].parent, input_tree->tree[i].children[0], input_tree->tree[i].children[1]);
//        // }
//
//        // create matrix cluster*leaves -- 0 if leaf is not in cluster, 1 if it is in cluster
//        long ** clusters = (long**)malloc((num_leaves - 1) * sizeof(long *));
//        for (long i = 0; i < num_leaves - 1; i++){
//            clusters[i] = (long*)malloc((num_leaves) * sizeof(long));
//        }
//
//        for (long i = 0; i < num_leaves ; i++){
//            for (long j = 0; j < num_leaves - 1; j++){
//                clusters[j][i] = 0; //initialise all entries to be 0
//            }
//            long j = i;
//            while (input_tree->tree[j].parent != -1){
//                j = input_tree->tree[j].parent;
//                clusters[j - num_leaves][i] = 1;
//            }
//            clusters[num_leaves - 2][i] = 1;
//        }
//
//        // convert matrix longo output string tree_str
//        sprintf(tree_str, "[{");
//        long tree_str_pos; //last position in tree_str that is filled with a character
//        for (long i = 0; i < num_leaves - 1; i++){
//            for (long j = 0; j < num_leaves; j++){
//                if (clusters[i][j] == 1){
//                    char leaf_str[num_digits_n + 1];
//                    sprintf(leaf_str, "%ld,", j+1);
//                    strcat(tree_str, leaf_str);
//                }
//            }
//            tree_str_pos = strlen(tree_str) - 1;
//            tree_str[tree_str_pos] = '\0'; // delete last comma
//            strcat(tree_str, "},{");
//            tree_str_pos +=2;
//        }
//        tree_str[tree_str_pos] = '\0'; // delete ,{ at end of tree_str
//        tree_str[tree_str_pos - 1] = '\0';
//        strcat(tree_str, "]");
//
//        for (long i = 0; i < num_leaves - 1; i++){
//            free(clusters[i]);
//        }
//        free(clusters);
//
//        return(tree_str);
//    }
//}

// FINDPATH without saving the path -- returns only the distance
//long findpath_distance(Tree *start_tree, Tree *dest_tree){
//    long num_leaves = start_tree->num_leaves;
//    long path_index = 0; // next position on path that we want to fill with a tree polonger
//    if (start_tree->tree == NULL){
//        printf("Error. Start tree doesn't exist.\n");
//    } else if (dest_tree->tree == NULL){
//        printf("Error. Destination tree doesn't exist.\n");
//    } else{
//        remove("./output/findpath.rtree");
//        // write_tree(start_tree->tree, num_leaves, "./output/findpath.rtree"); // this ruins the running time!!!!!!!!
//        long current_mrca; //rank of the mrca that needs to be moved down
//        Tree current_tree;
//        current_tree.tree = (Node*)malloc((2 * num_leaves - 1) * sizeof(Node));
//        current_tree.num_leaves = num_leaves;
//        for (long i = 0; i < 2 * num_leaves - 1; i++){
//            current_tree.tree[i] = start_tree->tree[i];
//        }
//        Tree * current_tree_polonger;
//        current_tree_polonger = &current_tree;
//        for (long i = num_leaves; i < 2 * num_leaves - 1; i++){
//            current_mrca = mrca(current_tree_polonger, dest_tree->tree[i].children[0], dest_tree->tree[i].children[1]);
//            // move current_mrca down
//            while(current_mrca != i){
//                bool did_nni = false;
//                for (long child_index = 0; child_index < 2; child_index++){ // find out if one of the children of current_tree.tree[current_mrca] has rank current_mrca - 1. If this is the case, we want to make an NNI
//                    if (did_nni == false && current_tree.tree[current_mrca].children[child_index] == current_mrca - 1){ // do nni if current longerval is an edge
//                        // check which of the children of current_tree.tree[current_mrca] should move up by the NNI move
//                        bool found_child = false; //indicate if we found the correct child
//                        long child_stays; // index of the child of current_tree.tree[current_mrca] that does not move up by an NNI move
//                        // find the index of the child of the parent of the node we currently consider -- this will be the index child_stays that we want in the end
//                        long current_child_index = dest_tree->tree[i].children[0]; // rank of already existing cluster in both current_tree.tree and dest_tree->tree
//                        while (found_child == false){
//                            while (current_tree.tree[current_child_index].parent < current_mrca - 1){ // find the x for which dest_tree->tree[i].children[x] is contained in the cluster induced by current_tree.tree[current_mrca - 1]
//                                current_child_index = current_tree.tree[current_child_index].parent;
//                            }
//                            // find the index child_stays
//                            if(current_tree.tree[current_child_index].parent == current_mrca - 1){
//                                found_child = true;
//                                if (current_tree.tree[current_tree.tree[current_child_index].parent].children[0] == current_child_index){
//                                    child_stays = 0;
//                                } else{
//                                    child_stays = 1;
//                                }
//                            } else{
//                                current_child_index = dest_tree->tree[i].children[1];
//                            }
//                        }
//                        nni_move(current_tree_polonger, current_mrca - 1, 1 - child_stays);
//                        did_nni = true;
//                        current_mrca--;
//                    }
//                }
//                if (did_nni == false){
//                    rank_move(current_tree_polonger, current_mrca - 1);
//                    current_mrca--;
//                }
//                path_index++;
//            }
//        }
//        free(current_tree.tree);
//    }
//    return path_index;
//}

//void test_random_single_edit(long n, long trials) {
//	long fail_count = 0;
//	for (long t = 0; t < trials; t++) {
//		Tree* T0 = (Tree*)malloc(sizeof(Tree));
//		generate_random_tree(n, T0);
//		Tree* T1 = (Tree*)malloc(sizeof(Tree));
//		copy_tree(T0, T1);
//		long* move = uniform_neighbour(T1);
//		if (move[0] == -1) {
//			continue;
//		}
//		Tree* R = (Tree*)malloc(sizeof(Tree));
//		generate_random_tree(n, R);
//
//		Precomp* precomp_data = precomp(T0, R);
//
//		Update* update_data = find_distance_diff(T0, R, move, precomp_data);
//		long diff = update_data->distance_diff;
//
//		long d0 = findpath_distance(T0, R);
//		long d1 = findpath_distance(T1, R);
//		if (diff != d1 - d0) {
//			fail_count += 1;
//		}
//	}
//	printf("failed %ld out of %ld\n", fail_count, trials);
//}

//void test_random_single_edit_chain(long n, long chain_length, long trials) {
//	long fail_total = 0;
//	for (long t = 0; t < trials; t++) {
//		long fail_count = 0;
//		Tree* T = (Tree*)malloc(sizeof(Tree));
//		generate_random_tree(n, T);
//		Tree* R = (Tree*)malloc(sizeof(Tree));
//		generate_random_tree(n, R);
//
//		Precomp* precomp_data = precomp(T, R);
//
//		Tree* T_current = (Tree*)malloc(sizeof(Tree));
//		copy_tree(T, T_current);
//		Tree* T_next = (Tree*)malloc(sizeof(Tree));
//		copy_tree(T_current, T_next);
//
//		for (long i = 0; i < chain_length; i++) {
//			long* move = uniform_neighbour(T_next);
//
//			Update* update_data = find_distance_diff(T_current, R, move, precomp_data);
//			long diff = update_data->distance_diff;
//
//			update_precomp(T_next, move, precomp_data, update_data);
//
//			copy_tree(T_next, T_current);
//		}
//		printf("failed %ld out of %ld\n", fail_count, chain_length);
//		if (fail_count > 0) fail_total++;
//	}
//	printf("total: failed %ld out of %ld\n", fail_total, trials);
//}

// ********** TESTING FUNCTIONS END **********

//long main() {
//	srand(time(NULL));
//	long n = 10;
//	long chain_length = 100;
//	long trials = 1000;
//	test_random_single_edit(n, trials);
//	test_random_single_edit_chain(n, chain_length, trials);
//}
