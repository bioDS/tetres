import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter


def plot_coords(coords, filename=None, colors=None, colorbar=False, centers=False, sizes=False, cmap=plt.cm.tab10, color_scaling=False, scale_centers=1):
    # Normalizing the colors

    internal_colors = colors
    if color_scaling:
        # this scales the colors to be in the interval 0 to 1
        # useful if you look at values such as likelihood or unnormed probabilities
        internal_colors = [(c - np.min(colors)) / (np.max(colors) - np.min(colors)) for c in colors]
    if colors is None:
        # no colors given, everything in the same color
        internal_colors = [1 for _ in range(coords.shape[0])]
    # the colors are
    internal_colors = [cmap(c) for c in internal_colors]

    if not sizes:
        sizes = [5 for _ in range(coords.shape[0])]

    if coords.shape[1] == 2:
        x, y = zip(*coords)
        plt.scatter(x, y, c=internal_colors, s=sizes, alpha=0.2)
        if centers:
            for i in range(1, len(centers) + 1):
                plt.scatter(x[-i], y[-i], color=internal_colors[-i], s=sizes[-i]*scale_centers, alpha=1)
                # adding a legend that identifies the centers with their respective colors
                markers = [plt.Line2D([0, 0], [0, 0], color=internal_colors[-i], marker='o', linestyle='') for i in
                           range(1, len(centers) + 1)]
                plt.legend(markers, centers, numpoints=1, loc="upper center", fancybox=True, ncol=len(centers),
                           bbox_to_anchor=(0.5, 1.15), shadow=True)
    elif coords.shape[1] == 3:
        x, y, z = zip(*coords)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, color=internal_colors, s=sizes, alpha=0.2)
        if centers:
            for i in range(1, len(centers) + 1):
                ax.scatter(x[-i], y[-i], z[-i], color=internal_colors[-i], s=sizes[-i]*scale_centers, alpha=1)
                # adding a legend that identifies the centers with their respective colors
                markers = [plt.Line2D([0, 0], [0, 0], color=internal_colors[-i], marker='o', linestyle='') for i in
                           range(1, len(centers) + 1)]
                plt.legend(markers, centers, numpoints=1, loc="upper center", fancybox=True, ncol=len(centers),
                           bbox_to_anchor=(0.5, 1.15), shadow=True)
    else:
        raise ValueError(f"Unsupported number of dimensions to plot {coords.shape[1]}!")

    if colorbar:
        # todo is this acurate?
        plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap=cmap), orientation="horizontal")

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename, dpi=400, format="png")
    plt.clf()
    plt.close()


def plot_density_over_coordinates(coords, density, filename, animate: bool = False):

    ###########################
    # TODO  This only makes sense for MDS coordinates because tsne is just to find clusters
    #       For the density it is important to have an implicit representaiton of the metric space not just a visual representation of probably close/far trees
    #
    ###########################
    if coords.shape[1] != 2:
        raise ValueError(f"Wrong dimenstion coordinates given {coords.shape[1]}!")
    if density.shape[0] != coords.shape[0]:
        raise ValueError(f"Wrong dimenstion density given, does not fit coords {density.shape[0]}\{coords.shape[0]}!")

    x, y = zip(*coords)
    z = density
    z = [(c - np.min(z)) / (np.max(z) - np.min(z)) for c in z]
    # todo there may be a problem with this scaling ?

    from scipy.interpolate import griddata

    # xi = np.linspace(np.min(x), np.max(x), 100)
    # yi = np.linspace(np.min(y), np.max(y), 100)

    xi, yi = np.mgrid[np.min(x):np.max(x):1000j, np.min(y):np.max(y):1000j]

    zi = griddata(points=(x, y), values=z, xi=(xi, yi), method='cubic', fill_value=np.min(z))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(x, y, z)
    dens = ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=plt.cm.Spectral_r)

    if animate:
        def init():
            ax.view_init(elev=15., azim=0)
            return [dens]
    
        def animate(i):
            ax.view_init(elev=15., azim=i * 45)
            return [dens]
    
        anim = FuncAnimation(fig, animate, init_func=init,
                             frames=8, interval=20000, blit=True)
        # # Save
        writergif = PillowWriter(fps=30)
        anim.save(f'{filename}/distribution.gif', writer=writergif)

    else:
        ax.view_init(elev=21., azim=45)
        plt.savefig(f"{filename}/distribution.png", dpi=200, bbox_inches="tight")
    
    plt.clf()
    return 0


# todo this should become part of the judgment package as this is not a visualization thing
def plot_all_chains_tsne(mchain, names=None, rf: bool = False):
    if names is None:
        names = [f"chain{i}" for i in range(mchain.m_MChains)]

    # labels = []

    cmap = plt.cm.Set1
    # todo add parameter to merge chains and see them as one ?
    # todo ore only plot specific chains like a list of indeces
    colors = []
    for i in range(mchain.m_MChains):
        colors.extend([cmap(i) for _ in range(len(mchain[i].trees))])
        # labels.append([names[i] for _ in range(len(mchain[i].trees))])

    combined_matrix = mchain.pwd_matrix_all(rf=rf)

        # print(combined_matrix.shape)
    # todo save the combined distance matrix!
    # todo make other visualization methods possible
    #  like metric MDS ! 
    # calling tsne for coordinate computaion
    # todo ideally this should be outsorced and possible with tsne recognizing the input differently maybe
    from treeoclock.visualize.tsne import _tsne_coords_from_pwd
    # todo it is totally fine to run tsne 10 times and select the solutin with the lowest KL divergence
    # todo the coords should be saved! and then loaded if called again!
    coords, kl_divergence = _tsne_coords_from_pwd(pwd_matrix=combined_matrix, dim=2)  # todo dimension

    # continue with the rest and color the dots accordingly

    # cmap = plt.cm.Spectral_r
    # colors = [cmap(c) for c in colors]

    # todo s is the size, so it would be possible to scale based on the colors
    #  or also another value
    s = [1 for _ in range(coords.shape[0])]

    if coords.shape[1] == 2:
        x, y = zip(*coords)
        plt.scatter(x, y, c=colors, s=s)
    elif coords.shape[1] == 3:
        x, y, z = zip(*coords)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c=colors, s=s)
    else:
        raise ValueError(f"Unsupported number of dimensions to plot {coords.shape[1]}!")

    # plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap=cmap), orientation="horizontal")
    # plt.legend([cmap(i) for i in ], names)

    markers = [plt.Line2D([0, 0], [0, 0], color=cmap(i), marker='o', linestyle='') for i in range(mchain.m_MChains)]
    plt.legend(markers, names, numpoints=1)
    
    plt.suptitle(f"KL-d: {kl_divergence}")
    
    # plt.show()
    plt.savefig(f"{mchain.working_dir}/plots/{mchain.name}_tsne_all_chains{'_rf' if rf else ''}.png", dpi=400, bbox_inches="tight")
    # plt.show()
    # if filename is None:
    #     plt.show()
    # else:
    #     plt.savefig(filename, dpi=400, format="png")
    plt.clf()
    plt.close()


# todo temporary funciton, to be merged with above which will allow different coordinates/ mds methods
def plot_multiple_chains_isomap(mchain, names):
    # names = [f"chain{i}" for i in range(coupled_chains.m_MChains)])  # todo this should be default names

    # labels = []
    combined_matrix = np.array([])
    cmap = plt.cm.Set1
    colors = []
    for i in range(mchain.m_MChains):
        colors.extend([cmap(i) for _ in range(len(mchain[i].trees))])
        # labels.append([names[i] for _ in range(len(mchain[i].trees))])
        cur_row = np.array([])
        for j in range(mchain.m_MChains):
            # print(i, j)
            if i < j:
                cur_row = np.concatenate((cur_row, mchain.pwd_matrix(i, j)),
                                         axis=1) if cur_row.size else mchain.pwd_matrix(i, j)
            elif i > j:
                cur_row = np.concatenate((cur_row, np.zeros(mchain.pwd_matrix(j, i).shape)),
                                         axis=1) if cur_row.size else np.zeros(mchain.pwd_matrix(j, i).shape)
            elif i == j:
                cur_row = np.concatenate((cur_row, mchain.pwd_matrix(i)),
                                         axis=1) if cur_row.size else mchain.pwd_matrix(i)
            # print(cur_row.shape)
        combined_matrix = np.concatenate((combined_matrix, cur_row)) if combined_matrix.size else cur_row
        # print(combined_matrix.shape)
    # todo save the combined distance matrix!
    # todo make other visualization methods possible
    #  like metric MDS !
    # todo ideally this should be outsorced and possible with tsne recognizing the input differently maybe

    from sklearn import manifold
    imap = manifold.Isomap(n_components=2, metric="precomputed")
    combined_matrix = combined_matrix + combined_matrix.transpose()
    coords = imap.fit_transform(combined_matrix)

    # continue with the rest and color the dots accordingly

    # cmap = plt.cm.Spectral_r
    # colors = [cmap(c) for c in colors]

    # todo s is the size, so it would be possible to scale based on the colors
    #  or also another value
    s = [5 for _ in range(coords.shape[0])]

    if coords.shape[1] == 2:
        x, y = zip(*coords)
        plt.scatter(x, y, c=colors, s=s)
    elif coords.shape[1] == 3:
        x, y, z = zip(*coords)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c=colors, s=s)
    else:
        raise ValueError(f"Not supported number of dimenstions to plot {coords.shape[1]}!")

    # plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap=cmap), orientation="horizontal")
    # plt.legend([cmap(i) for i in ], names)

    markers = [plt.Line2D([0, 0], [0, 0], color=cmap(i), marker='o', linestyle='') for i in range(mchain.m_MChains)]
    plt.legend(markers, names, numpoints=1)

    plt.savefig(f"{mchain.working_dir}/plots/{mchain.name}_isomap_all_chains.png", dpi=400, bbox_inches="tight")
    # plt.show()
    # if filename is None:
    #     plt.show()
    # else:
    #     plt.savefig(filename, dpi=400, format="png")
    plt.clf()
    plt.close()
