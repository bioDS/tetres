import networkx as nx
from scipy.sparse.csgraph import connected_components

import matplotlib.pyplot as plt
import numpy as np

# todo testing required!


def visualize_largest_cc(chain, i):
    dmi = chain.pwd_matrix(i)
    visualize_matrix(dmi, log_data=chain[i].log_data, outfile=f"{chain.working_dir}/{chain.name}_{i}_largest_cc.png")


def visualize_matrix(d_matrix, log_data, outfile):

    # todo problem if lcc_i is only one tree then this function has a problem, catch those cases and raise exception

    # c_dm = d_matrix.copy()
    # d_matrix[d_matrix > 1] = 0
    lcc_i = extract_largest_cc_indices(d_matrix.copy())
    c_dm = d_matrix[lcc_i,:][:,lcc_i]
    c_dm[c_dm > 1] = 0
    G = nx.from_numpy_matrix(c_dm)
    # G = nx.relabel_nodes(G, dict(zip(range(len(G.nodes())), range(len(G.nodes)))))
    # G = nx.drawing.nx_agraph.to_agraph(G)
    #
    # G.draw(outfile2, format='png', prog='neato')
    # G.draw(outfile, format='dot', prog='neato')


    # largest_cc = G.subgraph(max(nx.connected_components(G), key=len))
    largest_cc = G
    edge_list = [(u, v) for (u, v, d) in largest_cc.edges(data=True) if d['weight'] == 1]

    # todo add node size based on posterior probability

    # pos = nx.spring_layout(largest_cc, seed=7)
    pos = nx.nx_agraph.graphviz_layout(largest_cc, prog="fdp")

    # print(log_data.iloc[lcc_i])

    sizes = list(log_data.iloc[lcc_i]["likelihood"])
    min = np.min(sizes)
    max = np.max(sizes)
    # print(min, max, sizes)
    sizes = [(s - min)/(max - min) for s in sizes]

    nx.draw_networkx_nodes(largest_cc, pos, node_size=[s*100 for s in sizes])
    nx.draw_networkx_edges(largest_cc, pos, width=0.5, edgelist=edge_list)

    # nx.draw(G.subgraph(largest_cc))
    # plt.show()
    plt.savefig(fname=outfile,
                format="png", dpi=800)
    # .(outfile_png, format="png", prog="neato")

    # todo this largest cc or their indeces should also be somehow returned!
    # print("Starting")
    # nx.drawing.nx_agraph.write_dot(G.subgraph(largest_cc), outfile)
    # print("Finished!")


from collections import Counter


# todo testing required!
def extract_largest_cc_indices(d_matrix):
    d_matrix[d_matrix > 1] = 0
    ccs = connected_components(d_matrix)
    # print(len(ccs[1]))
    # print(len(set(ccs[1])))
    ccs_count = Counter(ccs[1])
    # print(ccs_count)
    val = max(ccs_count, key=lambda key: ccs_count[key])
    # print(val)
    ind = np.where(ccs[1] == val)[0]
    # print(ind)
    # print(len(ind))
    return ind

