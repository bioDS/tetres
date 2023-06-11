import matplotlib.pyplot as plt
import seaborn as sns


# todo the same for branch lengths somehow...
def get_lengths(tree):
    distances = []
    node2leaves = tree.etree.get_cached_content()
    index = 0
    for node in tree.etree.traverse('levelorder'):
        if len(node2leaves[node]) != 1:
            # if node.children:
            # if not node.is_leaf():
            if index == 0:
                node.name = node.dist
            else:
                # node.name = node.get_distance(tree_root))  # get_distance() is a very slow function
                node.name = node.dist + node.up.name  # Top down for loop makes this possible
            index += 1
            distances.append(node.name)
    distances = sorted(distances)

    bl = {}
    n = len(tree)
    for r in range(len(distances)-1):
        rank = n - 1 - r
        bl[rank] = distances[r+1] - distances[r]
    bl[1] = tree.etree.get_farthest_leaf()[1] - distances[-1]  # adding the rank 1 height separately,
    return bl
    # rank_etree = ctree_to_ete3(tree.ctree)
    # rc_dict = {}
    # for node in rank_etree.traverse("levelorder"):
    #     if len(node) > 1:
    #         c = node.children
    #         c0_leafs = []
    #         for leaf in c[0]:
    #             c0_leafs.append(int(leaf.name))
    #         c1_leafs = []
    #         for leaf in c[1]:
    #             c1_leafs.append(int(leaf.name))
    #         node_rank = int(node.name.replace("I", ""))
    #
    #         # c0 should always be the one with lowest integer among the two lists
    #         if min(c1_leafs) < min(c0_leafs):
    #             c1_leafs, c0_leafs = c0_leafs, c1_leafs
    #
    #         rc_dict[
    #             node_rank] = f"{','.join([str(i) for i in sorted(c0_leafs)])}|{','.join([str(i) for i in sorted(c1_leafs)])}"
    #
    # rc_dict[max(rc_dict.keys()) + 1] = ",".join(
    #     sorted([t for t in re.split(r"[,|]", rc_dict[max(rc_dict.keys())])], key=int))  # Setting root
    #
    #
    # for rank in sorted(rc_dict.keys(), reverse=True)[1:]:
    #     cur_split = rc_dict[rank]
    #     to_be_mod = rc_dict[rank + 1]
    #     to_be_swapped = ','.join(sorted([t for t in re.split(r'[,|]', cur_split)], key=int))
    #     rc_dict[rank] = re.sub(fr"{to_be_swapped}(?=$|\D)", cur_split, to_be_mod)
    # # return set of strings for a tree
    # return rc_dict




