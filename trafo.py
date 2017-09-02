import numpy as np
import networkx as nx
import random as rand

# allow degenerate coordinates
ALLOW_DEG_COORDS = True

#################
### testing   ###
#################


def plot_triangulation(G,
                       xmin=None,
                       xmax=None,
                       ymin=None,
                       ymax=None,
                       new_edges=None,
                       old_edges=None):
    import matplotlib.pyplot as plt
    from matplotlib import collections as mc

    xmin = G.graph['xmin'] if xmin is None else xmin
    xmax = G.graph['xmax'] if xmax is None else xmax
    ymin = G.graph['ymin'] if ymin is None else ymin
    ymax = G.graph['ymax'] if ymax is None else ymax

    # find edges that cross the domain due to periodic boundary conditions (heuristic)
    dx_noncross = (xmax - xmin) / 2.
    dy_noncross = (ymax - ymin) / 2.

    edges = []
    for e in G.edges():
        edges.append(G.get_edge_data(*e)['coords'])
    noncross_edges = [
        e for e in edges
        if abs(e[0][0] - e[1][0]) <= dx_noncross
        and abs(e[0][1] - e[1][1]) <= dy_noncross
    ]
    cross_edges = [e for e in edges if not e in noncross_edges]

    # make edges
    lc_noncross = mc.LineCollection(
        noncross_edges, colors='black', linewidths=2)
    lc_cross = mc.LineCollection(cross_edges, colors='0.8', linewidths=1)

    # plot edges
    ax = plt.axes()

    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))

    ax.add_collection(lc_cross)
    ax.add_collection(lc_noncross)

    # highlight?
    if new_edges:
        lines = [G[e[0]][e[1]]['coords'] for e in new_edges]
        lc = mc.LineCollection(lines, colors='red', linewidths=3)
        ax.add_collection(lc)

    if old_edges:
        lines = [(vertex_coord(G, e[0]), vertex_coord(G, e[1]))
                 for e in old_edges]
        lc = mc.LineCollection(
            lines, colors='blue', linestyle='dashed', linewidths=3)
        ax.add_collection(lc)

    plt.show()


def vertex_coord(G, n1):
    n2 = G.neighbors(n1)[0]
    n3 = G[n1][n2]['bordering'][0]
    return get_one_tricoord(G, n1, n2, n3)


#################
### checking  ###
#################


def get_fourcoords(graph, n1, n2, n3, n4):
    assert (n1, n3) in graph.edges() or (n3, n1) in graph.edges()
    c1, c2, c3 = get_tricoords(graph, n1, n2, n3)
    c1x, c4, c3x = get_tricoords(graph, n1, n4, n3)
    assert c1x == c1
    assert c3x == c3

    return c1, c2, c3, c4


def check_triangular(graph, n1, n2, n3):
    assert (n1 in graph.neighbors(n2)) == (n2 in graph.neighbors(n1))
    assert (n2 in graph.neighbors(n3)) == (n3 in graph.neighbors(n2))
    assert (n3 in graph.neighbors(n1)) == (n1 in graph.neighbors(n3))
    assert (n1 in graph.neighbors(n2)) and (n2 in graph.neighbors(n3)) and (
        n3 in graph.neighbors(n1))


def check_link_environment(graph, n1, n3):
    n2, n4 = graph[n1][n3]['bordering']

    # assert four different cells, i.e. no self-touching cells
    assert n1 != n2 != n3 != n4

    # assert that all cells have at least three neighbors
    assert len(graph.neighbors(n1)) > 2 and len(
        graph.neighbors(n2)) > 2 and len(graph.neighbors(n3)) > 2 and len(
            graph.neighbors(n4)) > 2

    # assert that no double triplets are present
    assert graph.degree(n1) + graph.degree(n2) > 6
    assert graph.degree(n2) + graph.degree(n3) > 6
    assert graph.degree(n3) + graph.degree(n4) > 6
    assert graph.degree(n4) + graph.degree(n1) > 6
    assert graph.degree(n1) + graph.degree(n3) > 6

    assert n1 in graph.neighbors(n2)
    assert n3 in graph.neighbors(n2)
    assert n1 in graph.neighbors(n4)
    assert n3 in graph.neighbors(n4)
    assert (n1, n2) in graph.edges() or (n2, n1) in graph.edges()
    assert (n2, n3) in graph.edges() or (n3, n2) in graph.edges()
    assert (n3, n4) in graph.edges() or (n4, n3) in graph.edges()
    assert (n4, n1) in graph.edges() or (n1, n4) in graph.edges()

    c1, c2, c3, c4 = get_fourcoords(graph, n1, n2, n3, n4)
    assert c1 in graph[n1][n3]['coords']
    assert c3 in graph[n1][n3]['coords']
    assert c2 in graph[n1][n2]['coords']
    assert c4 in graph[n1][n4]['coords']
    assert c2 in graph[n3][n2]['coords']
    assert c4 in graph[n3][n4]['coords']


#################
### interface ###
#################


def create_graph_from_vor(vor, periodic, xmin, xmax, ymin, ymax):
    """
    create graph from voronoi tesselation (periodic flag determines periodic
    boundary conditions; xmin/xmax/ymin/ymax are used for visualization only)
    """
    graph = nx.Graph()
    graph.graph['xmin'] = xmin
    graph.graph['xmax'] = xmax
    graph.graph['ymin'] = ymin
    graph.graph['ymax'] = ymax

    for n, v in enumerate(vor):
        neighbors = v['faces']
        if periodic: assert len(neighbors) > 2
        for neighbor in neighbors:
            adjacent = neighbor['adjacent_cell']
            bid1, bid2 = neighbor['vertices']
            b1 = [
                bn['adjacent_cell'] for bn in neighbors
                if bid1 in bn['vertices'] and not bid2 in bn['vertices']
            ]
            b2 = [
                bn['adjacent_cell'] for bn in neighbors
                if bid2 in bn['vertices'] and not bid1 in bn['vertices']
            ]
            assert len(b1) == len(b2) == 1
            b1 = b1[0]
            b2 = b2[0]
            coords = [vor[n]['original'], vor[adjacent]['original']]
            #border_coords = [vor[b1]['original'], vor[b2]['original']]
            if n >= 0 and adjacent >= 0 and b1 >= 0 and b2 >= 0:
                graph.add_edge(
                    n, adjacent, coords=coords,
                    bordering=[b1, b2])  #, border_coords=border_coords)
            else:
                assert not periodic

    if periodic:
        for e1, e2 in graph.edges():
            check_link_environment(graph, e1, e2)

    return graph


###############
### utility ###
###############


def four_edge(graph, e):
    edata = graph.get_edge_data(*e)
    n1, n3 = e
    n2, n4 = edata['bordering']
    return n1, n2, n3, n4


def get_one_tricoord(graph, n1, n2, n3):
    """return coordinates of vertex n1 assuming that (n1, n2, n3) form a (possibly degenerate) triangle."""
    c12 = graph[n1][n2]['coords']
    c13 = graph[n1][n3]['coords']
    c23 = graph[n2][n3]['coords']
    assert len(c12) == len(c13) == len(c23) == 2
    c = c12 + c13
    c.remove(c23[0])
    c.remove(c23[1])
    assert c[0] == c[1]
    return c[0]


def get_tricoords(graph, n1, n2, n3):
    """return list of coordinates of three vertices (n1, n2, n3) that form a triangle (in order)."""
    check_triangular(graph, n1, n2, n3)
    c1 = get_one_tricoord(graph, n1, n2, n3)
    c2 = get_one_tricoord(graph, n2, n3, n1)
    c3 = get_one_tricoord(graph, n3, n1, n2)
    if not ALLOW_DEG_COORDS:
        assert c1 != c2 != c3
    return [c1, c2, c3]


##################
### randomness ###
##################


def random_edge(graph):
    rindex = int(rand.uniform(0, len(graph.edges())))
    return graph.edges()[rindex]


def random_node(graph):
    rindex = int(rand.uniform(0, len(graph.nodes())))
    return graph.nodes()[rindex]


def random_multi(graph, periodic, largeness=7):
    '''pick random cell with at least largeness neighbors'''

    def condition(triplen, largeness):
        assert largeness is not None
        res = triplen < largeness
        return res

    # randomly select a multiple junction
    nnn = 0 if largeness is not None else 100  # number of neighbors
    count = 0
    while condition(nnn, largeness) and count < 1000:
        count += 1
        nc = random_node(graph)
        nn = graph.neighbors(nc)
        nnn = len(nn)

    if count == 1000:
        print 'lost in random_multi'
        ret = False
    else:
        ret = True
    return ret, nc


def random_triple(graph, nc):
    n1 = nc
    nn = graph.neighbors(nc)
    rand.shuffle(nn)
    n2 = nn[0]
    n3 = graph[n1][n2]['bordering'][0]
    return (n1, n2, n3)


#######################
### preferentiality ###
#######################


def preferential_flipartner(graph, nc, potentials):
    '''Choose a node x from potentials such that flip(nc, x) is the most preferential flip.'''

    wl = []
    dc = graph.degree(nc)
    for x in potentials:
        dx = graph.degree(x)
        b1, b2 = graph[x][nc]['bordering']
        d1, d2 = graph.degree(b1), graph.degree(b2)
        wl += [(d2 + d1, x)]

    nn = [y for (x, y) in sorted(wl)]
    return nn


def arbitrary_flipartner(graph, nc):
    '''Choose an arbitrary node forming an edge with nc.'''

    nn = graph.neighbors(nc)
    rand.shuffle(nn)
    return nn[0]


##############
### random ###
##############


def random_flip(graph, periodic):
    flip(graph, random_edge(graph), periodic)


def random_pref_flip(graph, periodic=True):
    nc = random_node(graph)
    bn = preferential_flipartner(graph, nc, graph.neighbors(nc))[0]
    flip(graph, (nc, bn), periodic)


def random_cellmerge(graph, preferential=True, periodic=True, debug=False):

    # get random node
    ret = True
    nc = random_node(graph)

    # get random neighboring cell as merge partner
    merge_partner = None

    # find the smallest neighbor, nd, of nc
    ndmin = min(graph.degree(x) for x in graph.neighbors(nc))
    nd = [x for x in graph.neighbors(nc) if graph.degree(x) == ndmin]
    assert len(nd) > 0
    merge_partner = nd[0]

    if ret:
        ret = cellmerge(
            graph,
            nc,
            preferential,
            periodic,
            merge_partners=[merge_partner],
            debug=debug)
        mpnew = max(graph.nodes()) + 1
        assert mpnew not in graph.nodes()
        atts = nx.get_edge_attributes(graph, 'bordering')
        for x, y in atts.keys():
            xx, yy = atts[(x, y)]
            if merge_partner == xx:
                graph[x][y]['bordering'] = [mpnew, yy]
            if merge_partner == yy:
                graph[x][y]['bordering'] = [xx, mpnew]
        nx.relabel_nodes(graph, mapping={merge_partner: mpnew}, copy=False)
    return ret


def random_celldivision(graph,
                        periodic=True,
                        preferential=True,
                        largeness=6,
                        debug=False):

    if largeness is None:
        ret = True
        nc = random_node(graph)
    else:
        ret, nc = random_multi(graph, periodic, largeness=largeness)

    if ret:
        newnode = random_triple(graph, nc)
        celldivision(
            graph, [nc], [newnode],
            periodic,
            preferential=preferential,
            debug=debug)
        ncnew = max(graph.nodes()) + 1
        assert ncnew not in graph.nodes()
        atts = nx.get_edge_attributes(graph, 'bordering')
        for x, y in atts.keys():
            xx, yy = atts[(x, y)]
            if nc == xx:
                graph[x][y]['bordering'] = [ncnew, yy]
            if nc == yy:
                graph[x][y]['bordering'] = [xx, ncnew]
        nx.relabel_nodes(graph, mapping={nc: ncnew}, copy=False)
    return ret


def random_doublecycle(graph,
                       periodic=True,
                       preferential=True,
                       collapsize=6,
                       debug=False):

    ret, nc = random_multi(graph, periodic, largeness=collapsize)
    if ret:
        ret = doublecycle(
            graph,
            nc,
            periodic=periodic,
            preferential=preferential,
            collapsize=collapsize,
            debug=debug)
    return ret


########################
### elemental trafos ###
########################


def tricollaps(graph, nc, periodic):
    '''collaps cell with 3 neighbors'''

    n1, n2, n3 = graph.neighbors(nc)
    check_link_environment(graph, n1, n3)
    check_triangular(graph, n1, n2, n3)

    c1, c2, cn1 = get_tricoords(graph, n1, n2, nc)
    c1, c3, cn2 = get_tricoords(graph, n1, n3, nc)
    assert cn1 == cn2
    cn = cn1
    assert get_tricoords(graph, n2, n3, nc) == [c2, c3, cn]

    # correct neighborhood relationships
    clist = {n1: c1, n2: c2, n3: c3}
    for x1, x2, y in [(n1, n2, n3), (n2, n3, n1), (n3, n1, n2)]:

        ## Switch neighbor identifier
        b0, b1 = graph[x1][x2]['bordering']
        if b0 == nc:
            graph[x1][x2]['bordering'][0] = y
        elif b1 == nc:
            graph[x1][x2]['bordering'][1] = y
        else:
            assert False

    # remove node (and automatically links)
    graph.remove_node(nc)
    check_link_environment(graph, n1, n2)
    check_link_environment(graph, n2, n3)
    check_link_environment(graph, n3, n1)
    check_triangular(graph, n1, n2, n3)


def tricreate(graph, triple=None, edge=None, periodic=True):

    if edge is not None:
        n1, n2, n3, n4 = four_edge(graph, edge)
    else:
        assert triple is not None
        n1, n2, n3 = triple

    check_link_environment(graph, n1, n3)
    check_triangular(graph, n1, n2, n3)

    c1, c2, c3 = get_tricoords(graph, n1, n2, n3)

    # create name for new node
    nnew = max(graph.nodes()) + 1
    assert not nnew in graph.nodes()

    # determine center of mass of triangle given by c1, c2, c3 as coordinates for new cell
    xnew = (c1[0] + c2[0] + c3[0]) / 3.
    ynew = (c1[1] + c2[1] + c3[1]) / 3.
    cnew = (xnew, ynew)
    if not ALLOW_DEG_COORDS:
        if cnew in [graph[x][y]['coords'][0]
                    for x, y in graph.edges()] or cnew in [
                        graph[x][y]['coords'][1] for x, y in graph.edges()
                    ]:
            cnew = (xnew + 0.1, ynew + 0.1)
        assert cnew not in [graph[x][y]['coords'][0] for x, y in graph.edges()]
        assert cnew not in [graph[x][y]['coords'][1] for x, y in graph.edges()]

    # correct neighborhood relationships
    cdict = {n1: c1, n2: c2, n3: c3}
    for x1, x2, y in [(n1, n2, n3), (n2, n3, n1), (n3, n1, n2)]:

        ## switch neighbor identifier
        b0 = graph[x1][x2]['bordering'][0]
        b1 = graph[x1][x2]['bordering'][1]
        if b0 == y:
            graph[x1][x2]['bordering'][0] = nnew
        elif b1 == y:
            graph[x1][x2]['bordering'][1] = nnew
        else:
            assert False

    # add triangular edges defining the new cell
    for (nrun, crun), (nb1, cb1), (nb2,
                                   cb2) in [((n1, c1), (n2, c2), (n3, c3)),
                                            ((n2, c2), (n3, c3), (n1, c1)),
                                            ((n3, c3), (n1, c1), (n2, c2))]:
        graph.add_edge(nnew, nrun, coords=[cnew, crun], bordering=[nb1, nb2])

    check_triangular(graph, n1, n2, n3)
    check_triangular(graph, n1, n2, nnew)
    check_triangular(graph, n1, n3, nnew)
    check_triangular(graph, n2, n3, nnew)
    check_link_environment(graph, n1, n2)
    check_link_environment(graph, n2, n3)
    check_link_environment(graph, n3, n1)
    check_link_environment(graph, n1, nnew)
    check_link_environment(graph, n2, nnew)
    check_link_environment(graph, n3, nnew)


def flip(graph, edge, periodic):
    # n1, n3 are original edge
    n1, n2, n3, n4 = four_edge(graph, edge)

    check_triangular(graph, n2, n1, n3)
    check_triangular(graph, n4, n1, n3)
    check_link_environment(graph, n1, n3)

    # make sure n1, n3 are not triangular to prevent a vanishing cell (in dual graph)
    triangle = len(graph.neighbors(n1)) <= 3 or len(graph.neighbors(n3)) <= 3

    # make sure no cells are created that touch more than once
    assert (not n4 in graph.neighbors(n2)) == (not n2 in graph.neighbors(n4))
    doubletouch = n4 in graph.neighbors(n2)

    if triangle or doubletouch:
        return False, None

    # determine coordinates of n2 and n4
    c2 = get_one_tricoord(graph, n2, n3, n1)
    c4 = get_one_tricoord(graph, n4, n1, n3)

    # flip edge
    graph.add_edge(n2, n4, coords=[c2, c4], bordering=[n1, n3])
    graph.remove_edge(n1, n3)

    # correct neighborhood relationships
    for x1, x2, y1, y2 in [(n1, n2, n3, n4), (n2, n3, n4, n1),
                           (n3, n4, n1, n2), (n4, n1, n2, n3)]:

        ## Switch neighbor identifier
        b0 = graph[x1][x2]['bordering'][0]
        b1 = graph[x1][x2]['bordering'][1]
        if b0 == y1:
            graph[x1][x2]['bordering'][0] = y2
        elif b0 == y2:
            graph[x1][x2]['bordering'][0] = y1
        elif b1 == y1:
            graph[x1][x2]['bordering'][1] = y2
        elif b1 == y2:
            graph[x1][x2]['bordering'][1] = y1
        else:
            assert False

    check_link_environment(graph, n2, n4)
    check_link_environment(graph, n1, n2)
    check_link_environment(graph, n2, n3)
    check_link_environment(graph, n3, n4)
    check_link_environment(graph, n4, n1)
    return True, (n2, n4)


###############################
### derived processes       ###
###############################


def multireduction(graph,
                   nc,
                   preferential,
                   periodic,
                   merge_partners=None,
                   spare_info=None,
                   debug=False):

    assert (spare_info is not None and merge_partners is not None) == False
    if spare_info is not None:
        spares, spare_target = spare_info
        assert len(spares) == 2
    nn = graph.neighbors(nc)

    if debug:
        print 'multireduction: entering'

    count = 0
    while len(nn) > 3 and count < 1000:
        count += 1
        graph_temp = graph

        if preferential:
            ret = False
            ntry = 0
            while not ret and ntry < 1:  #len(nn)/4.:

                if spare_info is not None:
                    pots = []
                    # ensure that newly created cells do not shrink by multireduction
                    if graph.degree(spares[0]) < spare_target[0]:
                        pots += [
                            x for x in nn
                            if spares[0] in graph[nc][x]['bordering']
                        ]
                    if graph.degree(spares[1]) < spare_target[1]:
                        pots += [
                            x for x in nn
                            if spares[1] in graph[nc][x]['bordering']
                        ]
                    # ensure that one of the newly created cells does not grow by multireduction
                    if graph.degree(
                            spares[1]) > spare_target[1] and spares[1] in nn:
                        pots += [spares[1]]
                    if pots == []:
                        pots = [x for x in nn if x not in spares]

                if merge_partners is not None:
                    pots = list(
                        set([
                            x
                            for sublist in [
                                graph[nc][mp]['bordering']
                                for mp in merge_partners if mp in nn
                            ] for x in sublist
                        ]))
                    if pots == []: pots = nn

                bn = preferential_flipartner(graph, nc, pots)[ntry]
                ret, new_edge = flip(graph, (nc, bn), periodic=periodic)
                ntry += 1
        else:
            rand.shuffle(nn)
            bn = nn[0]
            ret, new_edge = flip(graph, (nc, bn), periodic)

        if ret or not preferential:
            nn = graph.neighbors(nc)

        if len(nn) < 4:
            assert graph_temp == graph

        if debug:
            print 'multireduction: edge flipped'
            plot_triangulation(
                graph, new_edges=[new_edge], old_edges=[(nc, bn)])

    if count == 1000 and not preferential:
        print 'multireduction: lost in multireduction'
        return False
    elif count == 1000 and preferential:
        print 'multireduction: fixing'
        return multireduction(
            graph,
            nc,
            False,
            periodic,
            merge_partners=merge_partners,
            debug=debug)
    else:
        if debug:
            print 'multireduction: done'
        return True


def multicollaps(graph,
                 nc,
                 preferential,
                 periodic,
                 spare_info=None,
                 debug=False):
    multicollapsemerge(
        graph, nc, preferential, periodic, spare_info=spare_info, debug=debug)


def cellmerge(graph,
              nc,
              preferential,
              periodic,
              merge_partners=None,
              debug=False):
    multicollapsemerge(
        graph,
        nc,
        preferential,
        periodic,
        merge_partners=merge_partners,
        debug=debug)


def multicollapsemerge(graph,
                       nc,
                       preferential,
                       periodic,
                       merge_partners=None,
                       spare_info=None,
                       debug=False):

    if debug:
        print 'multicollapsemerge: initial graph'
        plot_triangulation(graph)

    multiok = multireduction(
        graph,
        nc,
        preferential,
        periodic,
        merge_partners=merge_partners,
        spare_info=spare_info,
        debug=debug)

    # remove node (and automatically links)
    if multiok:
        tricollaps(graph, nc, periodic=periodic)
        if debug:
            print 'multicollapsemerge: node removed'
            plot_triangulation(graph)

    if debug:
        print 'multicollapsemerge: done'

    return multiok


def celldivision(graph,
                 ncs,
                 newnodes,
                 periodic,
                 preferential=True,
                 debug=False):

    assert len(ncs) == len(newnodes)
    for k in range(len(ncs)):
        assert ncs[k] in newnodes[k] and ncs[0] in newnodes[k]

    n0 = ncs[0]
    nn = graph.neighbors(n0)
    nnn = len(nn)

    if debug:
        print 'celldivision: initial graph'
        plot_triangulation(graph)

    # create new node, perform checks and get indentifier of new node
    nnews = []
    for newnode in newnodes:
        nodes_temp = graph.nodes()
        tricreate(graph, triple=newnode, edge=None, periodic=periodic)
        nnew = [x for x in graph.nodes() if not x in nodes_temp]
        assert len(nnew) == 1
        nnew = nnew[0]
        assert (n0, nnew) in graph.edges() or (nnew, n0) in graph.edges()
        if debug:
            print 'celldivision: inserted new vertex'
            plot_triangulation(
                graph, new_edges=[(nnew, x) for x in graph.neighbors(nnew)])
        nnews += [nnew]

    # grow new node by edge flips flips until both cells have similar size
    ret = True
    for nnew, n01 in zip(nnews, ncs):
        nnn_temp = graph.degree(nnew)
        count = 0
        createinccond = np.floor(nnn / 2.) + 2
        while nnn_temp < createinccond and count < 1000:
            count += 1
            pots = [
                x for x in graph.neighbors(n01)
                if (nnew in graph[n01][x]['bordering'] and x not in nnews)
            ] if preferential else [graph[nnew][n01]['bordering'][0]]

            ret = False
            ntry = 0
            while not ret and ntry < 1:  #len(pots)/4.:
                bn = preferential_flipartner(graph, n01, pots)[ntry]
                ret, new_edge = flip(graph, (n01, bn), periodic=periodic)
                ntry += 1

            if ret:
                nnn_temp = graph.degree(nnew)
                if debug:
                    print 'celldivision: edge flipped'
                    plot_triangulation(
                        graph, new_edges=[new_edge], old_edges=[(n01, bn)])

        if count == 1000:
            print 'lost in celldivision'
            ret = False

    if debug:
        print 'celldivision: done'
    return [ret, nnews]


def doublecycle(graph,
                nc,
                periodic=True,
                preferential=True,
                collapsize=7,
                debug=False):
    '''creates two new nodes at both sides of link between nc and a six-sided neighbor ns,
     collapses nc and nl, the larger of the two cells bordering nc and ns
     interpretation: nc is an old neighboring cell, providing space for expasion, not the spawning cell'''

    ret = False

    # find a six sided neighboring cell: ns
    ns = [x for x in graph.neighbors(nc) if graph.degree(x) == 6]

    if len(ns) > 0:
        ns = ns[0]

        # find larger bordering cell: nl
        bordering = graph[nc][ns]['bordering']
        nl = bordering[0] if graph.degree(bordering[0]) > graph.degree(
            bordering[1]) else bordering[1]

        # get other bordering cell for complete description of the new nodes: nf
        nf1, nf2 = graph[nc][ns]['bordering']
        nf = nf1 if nf1 != nl else nf2

        # check the environments of the four cells nc, ns, nl, nf and of the triplets
        check_triangular(graph, nc, ns, nl)
        check_triangular(graph, nc, ns, nf)
        check_link_environment(graph, nc, nl)
        check_link_environment(graph, nc, ns)
        check_link_environment(graph, nc, nf)
        check_link_environment(graph, nf, ns)
        check_link_environment(graph, nl, ns)

        # if nl is large enough, perfom trafo
        if graph.degree(nl) > collapsize - 1:
            nnews = celldivision(
                graph, [nc, nl], [(nc, ns, nf), (nc, ns, nl)],
                periodic,
                preferential=preferential,
                debug=debug)[-1]
            Gintermed = graph.copy()
            target_sizes = [graph.degree(x) for x in nnews]
            multicollaps(
                graph,
                nc,
                preferential,
                periodic,
                spare_info=[nnews, target_sizes],
                debug=debug)
            multicollaps(
                graph,
                nl,
                preferential,
                periodic,
                spare_info=[nnews, target_sizes],
                debug=debug)
            ret = True
        else:
            Gintermed = graph.copy()
    else:
        Gintermed = graph.copy()
    return [ret, Gintermed]


def evolve_closed_cloud_cell(graph, periodic=False):
    random_celldivision(graph, periodic=periodic)
    random_cellmerge(graph, periodic=periodic)
    random_pref_flip(graph, periodic=periodic)
    random_pref_flip(graph, periodic=periodic)


def evolve_open_cloud_cell(graph, periodic=False):
    ret, Gintermed = random_doublecycle(graph, periodic=periodic)
    random_pref_flip(graph, periodic=periodic)
    random_pref_flip(graph, periodic=periodic)
    return Gintermed
