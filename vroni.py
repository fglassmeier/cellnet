import numpy as np
from scipy.stats import sem, linregress
import matplotlib.pyplot as plt
import matplotlib as mpl

font = {'weight': 'normal', 'size': 18}
mpl.rc('font', **font)
mpl.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble='\usepackage{sfmath}')


def get_degrees(G):
    return [G.degree(n) for n in G.nodes()]


def plot_degree_hist_from_deglist(deglist, axis=None, title=True):
    quantity = deglist
    func = np.abs(np.array(quantity) - 6.)**2
    sig = np.mean(func)
    bins = np.linspace(2.5, 9.5, 8)
    if axis is None:
        plt.hist(quantity, bins=bins, normed=True)
        plt.xlabel('number of neighbors')
        plt.ylabel('normalized frequency')
        if title: plt.title('$\sigma^2=$ %.2f' % sig)
        plt.show()
    else:
        axis.hist(quantity, bins=bins, normed=True)
        axis.set_xlabel('number of neighbors')
        axis.set_ylabel('normalized frequency')
        if title: axis.set_title('$\sigma^2=$ %.2f' % sig)


def nx_nneighbors_of_neighbors(G):
    return [(G.degree(k), np.average([G.degree(j) for j in G.neighbors(k)]))
            for k in G.nodes()]


def aboav_weaire(nneighbors_of_neighbors):
    xs = [d[0] for d in nneighbors_of_neighbors]
    ys = [d[0] * d[1] for d in nneighbors_of_neighbors]
    slope, intercept, r, p_value, std_err = linregress(xs, ys)
    return [xs, ys, slope, std_err, intercept, r, 6 - slope]


def plot_individual_aw(datlist, filename=None):
    xs, ys, slope, std_err, intercept, r, a = datlist
    cmap = plt.cm.gray_r
    plt.hist2d(
        xs,
        ys,
        cmap=cmap,
        bins=[[
            2.9, 3.1, 3.9, 4.1, 4.9, 5.1, 5.9, 6.1, 6.9, 7.1, 7.9, 8.1, 8.9,
            9.1, 9.9, 10.1, 100
        ],
              range(20, 100, 1)])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('number of data points')
    plt.plot(xs, [x * slope + intercept for x in xs])
    plt.title('$a=6-s=$ %.2f ($r=$ %.2f)' % (a, r))
    plt.xlabel('number of neighbors n')
    plt.ylabel('$n \cdot \langle N_n\\rangle$')
    x1, x2, y1, y2 = plt.axis()
    plt.axis((2.5, 10.5, 20, 60))
    if filename is not None:
        plt.savefig(filename, dpi=300)
    else:
        plt.show()
