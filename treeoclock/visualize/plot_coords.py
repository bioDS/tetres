
# todo in the future a cluster parameter will be added
# todo in the future likely a center parameter will be added to visualize summary trees
# todo in the future maybe adding a dimension parameter ?
import matplotlib.pyplot as plt
import numpy as np


def plot_coords(coords, filename=None, colors=None):

    # Normalizing the colors
    colors = [(c - np.min(colors)) / (np.max(colors) - np.min(colors)) for c in colors]
    cmap = plt.cm.Spectral_r
    colors = [cmap(c) for c in colors]

    # todo s is the size, so it would be possible to scale based on the colors
    #  or also another value
    s =[5 for _ in range(coords.shape[0])]

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

    plt.colorbar(plt.cm.ScalarMappable(norm=None, cmap=cmap), orientation="horizontal")

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename, dpi=400, format="png")
    plt.clf()


def plot_density_over_coordinates(coords, density, filename):

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


    zi = griddata(points=(x, y), values=z, xi=(xi, yi), method='linear', fill_value=np.min(z))


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # ax.scatter(x, y, z)
    ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=plt.cm.Spectral_r)
    # ax.plot_wireframe(xi, yi, zi, rstride=1, cstride=1, cmap="autumn")

    plt.savefig(filename, format="png", dpi=400)


    # ax.xaxis.pane.fill = False
    # ax.xaxis.pane.set_edgecolor('white')
    # ax.yaxis.pane.fill = False
    # ax.yaxis.pane.set_edgecolor('white')
    # ax.zaxis.pane.fill = False
    # ax.zaxis.pane.set_edgecolor('white')
    # ax.grid(False)
    # ax.w_zaxis.line.set_lw(0.)
    # ax.set_zticks([])
    #
    # start = 45
    # for i in range(8):
    #     ax.view_init(elev=15., azim=start*i)
    #     plt.savefig(f"{filename}/{i}.png", format="png", dpi=400)


    # plt.show()

    # from matplotlib import animation
    #
    # ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=plt.cm.Spectral_r)
    #
    # def init():
    #     # ax.scatter(xx, yy, zz, marker='o', s=20, c="goldenrod", alpha=0.6)
    #     ax.plot_surface(xi, yi, zi, rstride=1, cstride=1, cmap=plt.cm.Spectral_r)
    #     return fig,
    #
    # def animate(i):
    #     ax.view_init(elev=15., azim=i)
    #     return fig,
    #
    # # Animate
    # anim = animation.FuncAnimation(fig, animate, init_func=init, frames=360,
    #                                blit=True)
    # Save
    # plt.rcParams['animation.ffmpeg_path'] = '/home/lars/mambaforge/bin/ffmpeg'
    # anim.save(filename, writer="ffmpeg", extra_args=['-vcodec', 'libx264'])

    plt.clf()
    return 0


def plot_multiple_chains_tsne(mchain):

    # todo build new pw distance matrix from existing matrices
    combined_matrix = np.array([])
    colors = []
    for i in range(mchain.m_MChains):
        colors.extend([i for _ in range(len(mchain[i].trees))])
        cur_row = np.array([])
        for j in range(mchain.m_MChains):
            print(i, j)
            if i < j:
                cur_row = np.concatenate((cur_row, mchain.pwd_matrix(i, j)), axis=1) if cur_row.size else mchain.pwd_matrix(i, j)
            elif i > j:
                cur_row = np.concatenate((cur_row, np.zeros(mchain.pwd_matrix(j, i).shape)), axis=1) if cur_row.size else np.zeros(mchain.pwd_matrix(j, i).shape)
            elif i == j:
                cur_row = np.concatenate((cur_row, mchain.pwd_matrix(i)), axis=1) if cur_row.size else mchain.pwd_matrix(i)
            print(cur_row.shape)
        combined_matrix = np.concatenate((combined_matrix, cur_row)) if combined_matrix.size else cur_row
        # print(combined_matrix.shape)

    # todo save the combined distance matrix!

    # calling tsne for coordinate computaion
    # todo ideally this should be outsorced and possible with tsne recognizing the input differently maybe
    from treeoclock.visualize.tsne import _tsne_coords_from_pwd
    coords = _tsne_coords_from_pwd(pwd_matrix=combined_matrix, dim=2)  # todo dimension

    # continue with the rest and color the dots accordingly


    # cmap = plt.cm.Spectral_r
    # colors = [cmap(c) for c in colors]

    # todo s is the size, so it would be possible to scale based on the colors
    #  or also another value
    s =[5 for _ in range(coords.shape[0])]

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

    plt.show()
    # if filename is None:
    #     plt.show()
    # else:
    #     plt.savefig(filename, dpi=400, format="png")
    plt.clf()
