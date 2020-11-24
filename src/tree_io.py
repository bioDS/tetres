__author__ = 'Lena Collienne'
# Handling Tree input and output (to C)

import re
from call_findpath import *


# Read tree from string s in newick format -- assuming that the given tree is ultrametric!!
def read_newick(s):
    num_nodes = s.count(':') + 1 # number of nodes in input tree (internal nodes + leaves), root does not have ':'
    num_int_nodes = int((num_nodes - 1) / 2)
    num_leaves = int(num_nodes - num_int_nodes)

    tree_str = str(s) # copy input string -- we are going to manipulate it

    # While reading the newick string, we save the tree in the following way, converting it to C TREE afterwards
    int_node_index = 0 # index of current internal node -- NOT rank!
    parent_dict = dict() # for each node (leaves + internal) save index of parent (int_node_index)
    int_node_height = dict() # heights for all internal nodes (distance to leaves)

    while (len(tree_str) > 0): # recurse over tree and replace internal nodes with leaves labelled by internal_nodei for int i and save information about internal node height and children

        if int_node_index < num_int_nodes - 1: # as long as we don't reach the root
            pattern = r'\((\w+):(\[[^\]]*\])?((\d+.\d+E?\-?\d?)|(\d+)),(\w+):(\[[^\]]*\])?((\d+.\d+E?\-?\d?)|(\d+))\):(\[[^\]]*\])?((\d+.\d+E?\-?\d?)|(\d+))'
        else: # we reach the root -- string of form '(node1:x,node2:y)' left
            pattern = r'\((\w+):(\[[^\]]*\])?((\d+.\d+E?\-?\d?)|(\d+)),(\w+):(\[[^\]]*\])?((\d+.\d+E?\-?\d?)|(\d+))\);?'

        int_node_str = re.search(pattern, tree_str)

        # Save new internal node as parent of its two children
        parent_dict[int_node_str.group(1)] = int_node_index
        parent_dict[int_node_str.group(6)] = int_node_index

        int_node_height[int_node_index] = float(int_node_str.group(3)) # add height of new internal node to int_node_index   

        if int_node_index < num_int_nodes - 1: # insert new leaf 'internal_nodei' replacing the found pattern (pair of leaves)
            repl = "internal_node" + str(int_node_index) + ":" + str(float(int_node_height[int_node_index])+ float(int_node_str.group(12)))
            tree_str = re.sub(pattern, repl, tree_str, count = 1)
            int_node_index += 1
        else: # If we consider root, replace tree_str with empty str
            tree_str = ''

    # Now we use parent_list and int_node_height to create tree in C TREE data structure
    node_list = (NODE * num_nodes)()
    node_name_list = [] # Save the names of leaves in a list, such that entry i corresponds to leaf at position i in node list of tree structure

    # Initialise Node list
    for i in range(0, num_nodes):
        node_list[i].parent = -1
        node_list[i].children[0] = -1
        node_list[i].children[1] = -1
        node_name_list.append('')

    # sort internal nodes according to rank/times (increasing) -- output: list of keys (indices of internal nodes)
    int_node_list = [key for key, val in sorted(int_node_height.items(), key=lambda item:item[1])]

    leaf_index = 0 # Index for leaves, they are ordered alphabetical! (To make sure that two trees on same leaf set have same mapping of position in NODE array to label)
    parent_dict = dict(sorted(parent_dict.items(), key=lambda item: item[0])) # sort parent_dict alphabetical according to node names (internal_nodex + leaf labels) to make sure leaf labels are added to NODE list in alphabetical order

    for node in parent_dict:
        int_node_str = re.search(r'internal_node(\d*)', node)
        if int_node_str != None: # Consider internal node
            child_rank = num_leaves + int_node_list.index(int(int_node_str.group(1)))
            parent_rank = num_leaves + int_node_list.index(int(parent_dict[int_node_str.group()]))
            node_name_list[parent_rank] = str(node)
            # Set parent
            node_list[child_rank].parent = parent_rank
            # Set children (only one, make sure to use the slot children[0] ot children[1] that is not used already)
            if (node_list[parent_rank].children[0] == -1):
                node_list[parent_rank].children[0] = child_rank
            else:
                node_list[parent_rank].children[1] = child_rank
        else: # Consider leaf
            parent_rank = num_leaves + int_node_list.index(int(parent_dict[node]))
            # Set parent
            node_list[leaf_index].parent = parent_rank
            # Set children (only one, make sure to use the slot children[0] ot children[1] that is not used already)
            if (node_list[parent_rank].children[0] == -1):
                node_list[parent_rank].children[0] = leaf_index
            else:
                node_list[parent_rank].children[1] = leaf_index
            node_name_list[leaf_index] = str(node)
            leaf_index += 1

    output_tree = TREE(num_leaves, node_list)
    return(output_tree)

# Read trees from nexus file and save leaf labels as dict and trees as TREE_LIST
def read_nexus(file_handle):

    # To get the number of trees in the given file, we find the first and last line in the file that contain a tree (Line starts with tree)
    # Count number of lines in file
    number_of_lines = len(open(file_handle).readlines())
    last_line = number_of_lines

    # Find last line containing a tree
    for line in reversed(list(open(file_handle))):
        re_tree = re.search(r'^tree', line)
        if re_tree == None:
            last_line -= 1
        else: break

    # Find first line containing a tree
    first_line = 1
    for line in list(open(file_handle)):
        re_tree = re.search(r'^tree', line)
        if re_tree == None:
            first_line += 1
        else: break

    num_trees = last_line - first_line + 1 # Number of trees in nexus file

    # running variables for reading trees and displaying progress
    index = 0
    progress = 10

    name_dict = dict() # Save tree label names in dict
    trees = (TREE * num_trees)() # Save trees in an array to give to output TREE_LIST

    f = open(file_handle, 'r')

    leaf_labels = False
    # Save leaf labels -- returns empty dict if no abbreviation for leaf labels is used
    for line in f:
        if re.search(r'translate', line, re.I) != None:
            leaf_labels = True
        # Start reading leaf labels after 'translate' (signals start of this sequence)
        if leaf_labels == True:
            re_label = re.search(r'\s*(.+)\s(.+),', line)
            re_stop = re.search(r';', line)
            if re_stop != None:
                break
            elif re_label != None:
                name_dict[re_label.group(1)] = re_label.group(2)
    print(name_dict)

    # Read trees
    for line in f:
        re_tree = re.search(r'tree .* (\(.*\);)', line)
        if re_tree != None:
            current_tree = read_newick(re_tree.group(1))
            trees[index] = current_tree
            index += 1
            if int(100*index/num_trees) == progress:
                print(str(progress) + '% of trees are read')
                progress += 10

    tree_list = TREE_LIST(num_trees, trees)
    f.close()

    return(tree_list, name_dict)
