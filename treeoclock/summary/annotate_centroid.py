from treeoclock.trees.time_trees import TimeTreeSet, TimeTree
import numpy as np
import ete3


def annotate_centroid(tree: TimeTree, tree_set: TimeTreeSet):
    annotated_tree = tree
    n = annotated_tree.ctree.num_leaves
    rank_to_distances = {i: [] for i in range(n)}

    for t in tree_set:
        cur_distances = []
        rank_to_distances[0].append(t.etree.get_farthest_leaf()[1])
        node2leaves = t.etree.get_cached_content()
        index = 0
        for node in t.etree.traverse('levelorder'):
            if len(node2leaves[node]) != 1:
                if index == 0:
                    node.name = node.dist
                else:
                    node.name = node.dist + node.up.name  # Top down for loop makes this possible
                index += 1
                cur_distances.append(node.name)
        cur_distances = sorted(cur_distances)
        # sorted distances, cur_distances[1] is the height of node with rank n-2 i.e. the highest rank node apart from the root
        for i, el in enumerate(cur_distances):
            rank_to_distances[n-1-i].append(el)

    values = {k: np.mean(rank_to_distances[k]) for k in rank_to_distances.keys()}
    annotated_tree.etree = _ctree_to_ete3_annotation(annotated_tree.ctree, branch_lengths=values)
    return annotated_tree


def _ctree_to_ete3_annotation(ctree, branch_lengths):
    nl = ctree.num_leaves
    nn = (nl * 2) - 2  # number of nodes - 1, max index in ctree.tree

    def traverse(node):
        nonlocal branch_lengths
        nonlocal ctree
        nonlocal nl

        # Initializing current tree with the dist
        cur_t = ete3.Tree(dist=branch_lengths[ctree[node.children[1]].parent - nl + 1] - branch_lengths[node.parent - nl + 1],
                          name=f"I{ctree[node.children[1]].parent - nl + 1}")
        if node.children[0] >= nl:
            # Internal Node, setting dist as the difference between rank of parent and node itself
            cur_t.add_child(traverse(ctree[node.children[0]]))
        else:
            # Leaf node, dist is the rank of the parent node
            cur_t.add_child(name=node.children[0] + 1, dist=branch_lengths[0] - branch_lengths[ctree[node.children[0]].parent - nl + 1])
        if node.children[1] >= nl:
            # Internal Node, setting dist as the difference between rank of parent and node itself
            cur_t.add_child(traverse(ctree[node.children[1]]))
        else:
            # Leaf node, dist is the rank of the parent node
            cur_t.add_child(name=node.children[1] + 1, dist=branch_lengths[0] - branch_lengths[ctree[node.children[1]].parent - nl + 1])
        return cur_t

    t = traverse(ctree.tree[nn])
    return t
