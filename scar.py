import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt


def degree_difference_graph(Gnew, Gold):

    deg_new = Gnew.degree()
    deg_old = Gold.degree()
    DG = Gnew.copy()

    att = {}
    for x in Gnew.nodes():
        deg_diff = deg_new[x] - deg_old[x] if x in Gold.nodes() else deg_new[x]
        att[x] = deg_diff

    for x in Gold.nodes():
        if x not in Gnew.nodes():
            DG.add_node(x)
            att[x] = -deg_old[x]
            DG.add_edges_from(Gold.edges(x))

    assert len(att.keys()) == len(DG.nodes())

    nx.set_node_attributes(DG, 'degdiff', att)
    nx.set_node_attributes(DG, 'info', None)

    return DG


def degree_difference_stat(DG):

    sn = []
    ln = []
    for e1, e2 in DG.edges():
        ds = DG.node[e1]['degdiff']
        dl = DG.node[e2]['degdiff']
        if dl < ds: ds, dl = dl, ds
        if ds != 0 or dl != 0:
            sn += [ds]
            ln += [dl]

    return (sn, ln)


def plot_degree_difference_histogram(smallneigh,
                                     largeneigh,
                                     filename=None,
                                     ax=None,
                                     log=False):
    if ax is not None: assert filename is None
    norm = mpl.colors.LogNorm() if log else None
    vmin = 1e-3 if log else 0
    vmax = 1e-1 if log else 0.05
    bins = np.arange(-9.5, 10., 1.)
    if ax is not None:
        im = ax.hist2d(
            smallneigh,
            largeneigh,
            bins=[bins, bins],
            normed=True,
            vmin=vmin,
            vmax=vmax,
            norm=norm)
        ax.set_xlabel('$\Delta n_1$')
        ax.set_ylabel('$\Delta n_2$')
    else:
        im = plt.hist2d(
            smallneigh,
            largeneigh,
            bins=[bins, bins],
            normed=True,
            vmin=vmin,
            vmax=vmax,
            norm=norm)
        plt.xlabel('$\Delta n_1$')
        plt.ylabel('$\Delta n_2$')
    plt.colorbar(im[-1])
    if filename is not None:
        plt.savefig(filename)
    elif ax is None:
        plt.show()
