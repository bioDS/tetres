import networkx as nx

# todo testing required!
def visualize_matrix(d_matrix, outfile, outfile_png):
    d_matrix[d_matrix > 1] = 0
    G = nx.from_numpy_matrix(d_matrix)
    # G = nx.relabel_nodes(G, dict(zip(range(len(G.nodes())), range(len(G.nodes)))))
    # G = nx.drawing.nx_agraph.to_agraph(G)
    #
    # G.draw(outfile2, format='png', prog='neato')
    # G.draw(outfile, format='dot', prog='neato')

    largest_cc = max(nx.connected_components(G), key=len)
    # G.subgraph(largest_cc).(outfile_png, format="png", prog="neato")

    print("Starting")
    nx.drawing.nx_agraph.write_dot(G.subgraph(largest_cc), outfile)
    print("Finished!")


