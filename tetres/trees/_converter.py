import sys
import ete3
from tetres.trees._ctrees import NODE, TREE

# TODO Documentation


def ete3_to_ctree(tree):
    distances = []
    node2leaves = tree.get_cached_content()
    # tree_root = tree.get_tree_root()

    index = 0
    for node in tree.traverse('levelorder'):
        if len(node2leaves[node]) != 1:
            # if node.children:
            # if not node.is_leaf():
            if index == 0:
                node.name = node.dist
            else:
                # node.name = node.get_distance(tree_root)  # get_distance() is a very slow function
                node.name = node.dist + node.up.name  # Top down for loop makes this possible
            index += 1
            distances.append(node.name)

    # TODO node.time = ? needs to be set to the correct rank of each node!
    num_nodes = len(node2leaves)  # Number of nodes
    num_leaves = int(((num_nodes - 1) / 2) + 1)
    node_list = (NODE * num_nodes)()

    distances = sorted(distances)
    # if not len(distances.keys()) == num_leaves - 1:
    if not len(set(distances)) == num_leaves - 1:
        # TODO argument so that the ranked trees will be ignored fully ...
        sys.exit('Distances to root not unique! \n'
                 'This has to be resolved!')
    for node in tree.traverse('levelorder'):
        if len(node2leaves[node]) == 1:
            # if not node.children:
            # if node.is_leaf():
            node_list[int(node.name) - 1].parent = num_nodes - (distances.index(node.up.name) + 1)
        else:
            if node.name == 0.0:
                node_list[num_nodes - (distances.index(node.name) + 1)].parent = num_nodes - 1
            else:
                node_list[num_nodes - (distances.index(node.name) + 1)].parent = \
                    num_nodes - (distances.index(node.up.name) + 1)
            # current_children = node.get_children()  # get_children() is slow
            current_children = node.children
            if len(current_children) != 2:
                sys.exit('Not a binary tree!')

            # Child 0
            if len(node2leaves[current_children[0]]) == 1:
                # if not current_children[0].children:
                # if current_children[0].is_leaf():
                node_list[num_nodes - (distances.index(node.name) + 1)].children[0] = \
                    int(current_children[0].name) - 1
            else:
                node_list[num_nodes - (distances.index(node.name) + 1)].children[0] = \
                    num_nodes - (distances.index(current_children[0].name) + 1)
            # Child 1
            if len(node2leaves[current_children[1]]) == 1:
                # if not current_children[1].children:
                # if current_children[1].is_leaf():
                node_list[num_nodes - (distances.index(node.name) + 1)].children[1] = \
                    int(current_children[1].name) - 1
            else:
                node_list[num_nodes - (distances.index(node.name) + 1)].children[1] = \
                    num_nodes - (distances.index(current_children[1].name) + 1)
    return TREE(num_leaves, node_list, -1)


def ctree_to_ete3(ctree):
    nl = ctree.num_leaves
    nn = (nl * 2) - 2  # number of nodes - 1, max index in ctree.tree

    def traverse(node):
        nonlocal ctree
        nonlocal nl

        # Initializing current tree with the dist
        cur_t = ete3.Tree(dist=node.parent - ctree[node.children[1]].parent,
                          name=f"I{ctree[node.children[1]].parent - nl + 1}")
        if node.children[0] >= nl:
            # Internal Node, setting dist as the difference between rank of parent and node itself
            cur_t.add_child(traverse(ctree[node.children[0]]))
        else:
            # Leaf node, dist is the rank of the parent node
            cur_t.add_child(name=node.children[0] + 1, dist=ctree[node.children[0]].parent-nl+1)
        if node.children[1] >= nl:
            # Internal Node, setting dist as the difference between rank of parent and node itself
            cur_t.add_child(traverse(ctree[node.children[1]]))
        else:
            # Leaf node, dist is the rank of the parent node
            cur_t.add_child(name=node.children[1] + 1, dist=ctree[node.children[1]].parent-nl+1)
        return cur_t
    t = traverse(ctree.tree[nn])
    return t
