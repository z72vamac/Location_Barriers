# TSPN-B

import time
import itertools

import gurobipy as gp
from gurobipy import GRB
from matplotlib.patches import Circle

import auxiliar_functions as af
from data import *
import neighborhood as neigh


def h_kmedian_n(barriers, sources, targets, k, wL = 50, prepro=True, log=False, picture=False, time_limit=7200, init=False):

    S = len(sources)
    T = len(targets)

    first_time = time.time()

    # Indices of the vertices identified with the sources
    vertices_source = list(itertools.product(range(1, len(sources) + 1), range(1)))

    # Indices of the vertices identified with the vertices of the barriers
    vertices_barrier = list(itertools.product(range(1000, 1000 + len(barriers)), range(2)))

    # Indices of the vertices identified with the sources
    vertices_target = list(itertools.product(range(-len(targets), 0), range(1)))
    vertices_target = vertices_target[::-1]


    # Indices of the edges joining a source and a barrier
    edges_source = []

    for (a, b) in vertices_source:
        for c, d in vertices_barrier:
            if prepro:
                # Point of the barrier to check if is visible by the neighborhood
                point = barriers[c-1000][d]

                # Neighborhood to check if it is visible by the point
                neighborhood = sources[a - 1]

                if af.cansee(point, neighborhood, barriers):
                    # Appending the feasible edges to edges_neighborhood
                    edges_source.append((a, b, c, d))
                    # edges_neighborhood.append((c, d, a, b))
            else:
                edges_source.append((a, b, c, d))
                # edges_neighborhood.append((c, d, a, b))

    edges_barrier = []
    for v, i in vertices_barrier:
        for w, j in vertices_barrier:
            if v != w:
                if prepro:
                    barrier = [barriers[v-1000][i], barriers[w-1000][j]]

                    if any([af.intersect(barrieri, barrier) for barrieri in barriers]):
                        pass
                    else:
                        edges_barrier.append((v, i, w, j))
                else:
                    edges_barrier.append((v, i, w, j))

    indices_barriers = [(v, 0, v, 1) for v in range(1000, 1000 + len(barriers))]

    edges_target = []

    for (a, b) in vertices_target:
        for c, d in vertices_barrier:
            if prepro:
                # Point of the barrier to check if is visible by the neighborhood
                point = barriers[c-1000][d]

                # Neighborhood to check if it is visible by the point
                neighborhood = targets[abs(a) - 1]

                if af.cansee(point, neighborhood, barriers):
                    # Appending the feasible edges to edges_neighborhood
                    # edges_neighborhood.append((a, b, c, d))
                    edges_target.append((c, d, a, b))
            else:
                # edges_neighborhood.append((a, b, c, d))
                edges_target.append((c, d, a, b))


    vertices_total = vertices_source + vertices_barrier + vertices_target
    edges_total = edges_source + edges_barrier + edges_target

    if log:
        print("vertices_source = " + str(vertices_source))
        print("vertices_barrier = " + str(vertices_barrier))
        print("vertices_target = " + str(vertices_target))


        print("edges_source = " + str(edges_source))
        print("edges_barrier = " + str(edges_barrier))
        print("edges_free = " + str(edges_target))

    epsilon_index = []  # epsilon(S / T, B, i) = 1 si (P_{S/T}, P_B^i)\in E_{S/T}

    # print(epsilon_index)

    point_index = []
    for a, b in vertices_source+vertices_target:
        for dim in range(2):
            point_index.append((a, b, dim))


    # for a, b, c, d in edges_total:
    #     for e, f in vertices_neighborhood:
    #         y_index.append((a, b, c, d, e))

    alpha_index = []

    # for a, b in vertices_total:
    #     for c, d in vertices_total:
    #         for e, f in vertices_total:
    #             if (c, d, e, f) in indices_barriers or a > 0:
    #                 alpha_index.append((a, b, c, d, e, f))

    for a, b in vertices_barrier:
        for c, d, e, f in edges_total:
            alpha_index.append((a, b, c, d, e, f))

    for a, b in vertices_total:
        for (c, d, e, f) in indices_barriers:
            alpha_index.append((a, b, c, d, e, f))

    if log:
        print("alpha = " + str(alpha_index))

    beta_index = []

    # for a, b in vertices_total:
    #     for c, d in vertices_total:
    #         for e, f in vertices_total:
    #             for g, h in vertices_total:
    #                 if ((a, b, c, d) in indices_barriers and (e, f, g, h) in edges_total) or ((a, b, c, d) in edges_total and (e, f, g, h) in indices_barriers):
    #                     beta_index.append((a, b, c, d, e, f, g, h))

    for a, b, c, d in indices_barriers:
        for e, f, g, h in edges_total:
            beta_index.append((a, b, c, d, e, f, g, h))
            beta_index.append((e, f, g, h, a, b, c, d))

    # print([index for index in beta_index if (index in list(set(beta_index)))])
    # beta_index = list(beta_index)-list(set(beta_index))


    if log:
        print("beta = " + str(beta_index))

    gamma_index = beta_index

    delta_index = gamma_index

    epsilon_index = []

    for a, b, c, d in edges_total:
        epsilon_index.append((a, b, c, d))


    x_index = edges_total

    aux_index = []

    for a, b, c, d in edges_total:
        for e, f in vertices_target:
            aux_index.append((a, b, c, d, e, f))


    # for e, f in vertices_source:
    #     for g, h in vertices_target:
    #         for a, b, c, d in edges_total:
    #             if a == e or c == g or (a, b, c, d) in edges_barrier:
    #                 x_index.append((a, b, c, d, e, f, g, h))

    # print(x_index)

    y_index = vertices_source


    z_index = []

    for a, b in vertices_source:
        for c, d in vertices_target:
            z_index.append((a, b, c, d))

    if log:
        print("x_index = " + str(x_index))
        print("y_index = " + str(y_index))
        print("z_index = " + str(z_index))

    # P_S and P_T: indices of the points in the neighborhoods
    # p_index = []
    # for a, b, c, d in edges_total:
    #     for e, f in vertices_neighborhood:
    #         p_index.append((a, b, c, d, e))
    p_index = x_index

    paux_index = aux_index

    # for index in vertices_neighborhood:
    #     for dim in range(2):
    #         p_index.append((index[0], index[1], dim))

    if log:
        print("p_index = " + str(p_index))

    dist_index = edges_total

    if log:
        print("dist = " + str(dist_index))

    dif_index = []

    for a, b, c, d in edges_total:
        for dim in range(2):
            dif_index.append((a, b, c, d, dim))

    if log:
        print("dif = " + str(dif_index))

    # socp variables:
    d_inside_index = vertices_source + vertices_target

    dif_inside_index = []

    for (a, b) in vertices_source + vertices_target:
        for dim in range(2):
            dif_inside_index.append((a, b, dim))

    if log:
        print("d_inside_index = " + str(d_inside_index))
        print("dif_inside_index = " + str(dif_inside_index))

    # z variables:
    # z_index = list(itertools.product(vertices_neighborhood, vertices_neighborhood))
    # print("z_index = " + str(z_index))

    # f variables:
    # g_index = y_index
    #
    # if log:
    #     print("g_index = " + str(g_index))

    model = gp.Model('Model: H-KMedian-N')

    epsilon = model.addVars(epsilon_index, vtype=GRB.BINARY, name='epsilon')
    p = model.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name='p')
    paux = model.addVars(paux_index, vtype=GRB.CONTINUOUS, lb=0.0, name = 'paux')
    x = model.addVars(x_index, vtype=GRB.INTEGER, lb=0.0, ub = len(targets), name = 'x')
    aux = model.addVars(aux_index, vtype = GRB.BINARY, name = 'aux')
    y = model.addVars(y_index, vtype=GRB.BINARY, name='y')
    z = model.addVars(z_index, vtype=GRB.BINARY, name = 'z')
    dist = model.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dist')
    dif = model.addVars(dif_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')
    delta = model.addVars(delta_index, vtype=GRB.BINARY, name='delta')
    gamma = model.addVars(gamma_index, vtype=GRB.BINARY, name='gamma')
    # beta = model.addVars(beta_index, vtype = GRB.BINARY, name = 'beta')
    # alpha = model.addVars(alpha_index, vtype = GRB.BINARY, name = 'alpha')

    alpha = model.addVars(alpha_index, vtype=GRB.BINARY, name='alpha')

    beta = model.addVars(beta_index, vtype=GRB.BINARY, name='beta')

    # point = model.addVars(p_index, vtype = GRB.CONTINUOUS, lb = 0.1, ub = 99.9, name = 'point')
    point = model.addVars(point_index, vtype=GRB.CONTINUOUS, name='point')
    # point_star = model.addVars(point_star_index, vtype=GRB.CONTINUOUS, name='point_star')

    d_inside = model.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='d_inside')
    dif_inside = model.addVars(dif_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif_inside')
    # landa = model.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='landa')
    # z = model.addVars(z_index, vtype = GRB.BINARY, name = 'z')
    # g = model.addVars(g_index, vtype=GRB.CONTINUOUS, lb=0.0, name='g')

    model.update()

    if init:
        time_h, objval_h = heuristic(barriers, neighborhoods)

    # Alpha constraint
    L = -100000
    U = 100000

    # alpha-C
    for a, b, c, d, e, f in alpha_index:
        # Dos primeras
        # print((a, b, c, d, e, f))
        L = -100000
        U = 100000
        if (c, d, e, f) in indices_barriers:
            if (a, b) in vertices_source:
                # L, U = af.estima_det(sources[a - 1], [barriers[c - 1000][0], barriers[c - 1000][1], barriers[e - 1000][0], barriers[e - 1000][1]])
                model.addConstr((1-alpha[a,b,c,d,e,f])*L <= af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c-1000][d], barriers[e-1000][f]))
                model.addConstr(-U*alpha[a,b,c,d,e,f]<= -af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c-1000][d], barriers[e-1000][f]))
            elif (a, b) in vertices_barrier:
                model.addConstr((1-alpha[a,b,c,d,e,f])*L <= af.determinant(barriers[a-1000][b], barriers[c-1000][d], barriers[e-1000][f]))
                model.addConstr(-U*alpha[a,b,c,d,e,f]<= -af.determinant(barriers[a-1000][b], barriers[c-1000][d], barriers[e-1000][f]))
            else:
                # L, U = af.estima_det(targets[abs(a) - 1], [barriers[c - 1000][0], barriers[c - 1000][1], barriers[e - 1000][0], barriers[e - 1000][1]])
                model.addConstr((1-alpha[a,b,c,d,e,f])*L <= af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c-1000][d], barriers[e-1000][f]))
                model.addConstr(-U*alpha[a,b,c,d,e,f]<= -af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c-1000][d], barriers[e-1000][f]))

        else:
            if (c, d, e, f) in edges_source:
                model.addConstr((1-alpha[a,b,c,d,e,f])*L <= af.determinant(barriers[a-1000][b], [point[c, d, 0], point[c, d, 1]], barriers[e-1000][f]))
                model.addConstr(af.determinant(barriers[a-1000][b], [point[c, d, 0], point[c, d, 1]], barriers[e-1000][f]) <= U*alpha[a,b,c,d,e,f])
            elif (c, d, e, f) in edges_barrier:
                model.addConstr((1-alpha[a,b,c,d,e,f])*L <= af.determinant(barriers[a-1000][b],  barriers[c-1000][d], barriers[e-1000][f]))
                model.addConstr(U*alpha[a,b,c,d,e,f] >= af.determinant(barriers[a-1000][b],  barriers[c-1000][d], barriers[e-1000][f]))
            else:
                model.addConstr((1-alpha[a,b,c,d,e,f])*L <= af.determinant(barriers[a-1000][b], barriers[c-1000][d], [point[e, f, 0], point[e, f, 1]]))
                model.addConstr(af.determinant(barriers[a-1000][b], barriers[c-1000][d], [point[e, f, 0], point[e, f, 1]]) <= U*alpha[a,b,c,d,e,f])

    model.update()

    for a, b, c, d, e, f, g, h in beta_index:
        if ((a, b, c, d) in indices_barriers) or ((e, f, g, h) in indices_barriers):
            model.addConstr(beta[a, b, c, d, e, f, g, h] == 2*gamma[a, b, c, d, e, f, g, h] - alpha[a, b, e, f, g, h] - alpha[c, d, e, f, g, h] + 1)

            model.addConstr(gamma[a, b, c, d, e, f, g, h] <= alpha[a, b, e, f, g, h])
            model.addConstr(gamma[a, b, c, d, e, f, g, h] <= alpha[c, d, e, f, g, h])
            model.addConstr(gamma[a, b, c, d, e, f, g, h] >= alpha[a, b, e, f, g, h] + alpha[c, d, e, f, g, h] - 1)

        if (e, f, g, h) in indices_barriers:
            model.addConstr(0.5*(beta[a, b, c, d, e, f, g, h] + beta[e, f, g, h, a, b, c, d]) <= delta[a, b, c, d, e, f, g, h])
            model.addConstr(2 * (beta[a, b, c, d, e, f, g, h] + beta[e, f, g, h, a, b, c, d]) >= delta[a, b, c, d, e, f, g, h])

    for a, b, c, d in epsilon_index:
        model.addConstr(gp.quicksum(delta[a1, b1, c1, d1, e, f, g, h] for a1, b1, c1, d1, e, f, g, h in delta_index if (a1==a and b1==b and c1==c and d1==d and (e, f, g, h) in indices_barriers)) - len(barriers) + 1 <= epsilon[a, b, c, d])
        model.addConstr(len(barriers) * epsilon[a, b, c, d] <= gp.quicksum(delta[a1, b1, c1, d1, e, f, g, h] for a1, b1, c1, d1, e, f, g, h in delta_index if (a1==a and b1==b and c1==c and d1==d and (e, f, g, h) in indices_barriers)))

        model.addConstr(x[a, b, c, d] <= T*epsilon[a, b, c, d])

    model.update()

    # epsilon-C and y-C
    # for a, b, c, d in epsilon.keys():
    #     # print((a, b, c, d))
    #     # model.addConstr((delta.sum(a, b, c, d, '*') - len(barriers)) + 1 <= epsilon[a, b, c, d])
    #     # model.addConstr(len(barriers) * epsilon[a, b, c, d] <= delta.sum(a, b, c, d, '*'))
    #
    #     # if a < 0 and b >= 0 and c >= 0:
    #     for e, f in vertices_neighborhood:
    #         model.addConstr(y[a, b, c, d, e] + y[c, d, a, b, e] <= 2 * epsilon[a, b, c, d])

        # if a < 0 and b < 0 and c == 0:
        #     model.addConstr(y[a, b, c] + y[b, a, c] <= 2*epsilon[a, b, c])

        # if no_ve(neighborhoods[abs(a)-1], barriers[b][0]) and no_ve(neighborhoods[abs(a)-1], barriers[b][1]):
        #     model.addConstr(epsilon[a, b, c] <= 0)

    # NS and NT constraints
    for a, b, dim in point.keys():

        if (a, b) in vertices_source:
            source = sources[a - 1]

            # model.addConstr(point[a, b, dim] == source.center[dim])
            if type(source) is neigh.Circle:
                model.addConstr(dif_inside[a, b, dim] >= point[a, b, dim] - source.center[dim])
                model.addConstr(dif_inside[a, b, dim] >= source.center[dim] - point[a, b, dim])

                model.addConstr(gp.quicksum(dif_inside[a, b, dim] * dif_inside[a, b, dim] for dim in range(2)) <= d_inside[a, b] * d_inside[a, b])
                model.addConstr(d_inside[a, b] <= source.radii)

            if type(source) is neigh.Poligonal:
                model.addConstrs(point[a, b, dim] == landa[a, b] * source.V[0][dim] + (1 - landa[a, b]) * source.V[1][dim] for dim in range(2))

        elif (a, b) in vertices_target:
            target = targets[abs(a) - 1]

            # model.addConstr(point[a, b, dim] == target.center[dim])

            if type(target) is neigh.Circle:
                model.addConstr(dif_inside[a, b, dim] >= point[a, b, dim] - target.center[dim])
                model.addConstr(dif_inside[a, b, dim] >= target.center[dim] - point[a, b, dim])

                model.addConstr(gp.quicksum(dif_inside[a, b, dim] * dif_inside[a, b, dim] for dim in range(2)) <= d_inside[a, b] * d_inside[a, b])
                model.addConstr(d_inside[a, b] <= target.radii)

            if type(target) is neigh.Poligonal:
                model.addConstrs(point[a, b, dim] == landa[a, b] * target.V[0][dim] + (1 - landa[a, b]) * target.V[1][dim] for dim in range(2))

    # dist constraints
    for a, b, c, d, dim in dif_index:

        if (a, b, c, d) in edges_barrier:
            model.addConstr(dist[a, b, c, d] == np.linalg.norm(np.array(barriers[a-1000][b]) - np.array(barriers[c-1000][d])))

        if (a, b, c, d) in edges_source:
            model.addConstr(dif[a, b, c, d, dim] >= point[a, b, dim] - barriers[c-1000][d][dim])
            model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + barriers[c-1000][d][dim])
            model.addConstr(gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] * dist[a, b, c, d])

        elif (a, b, c, d) in edges_target:
            model.addConstr(dif[a, b, c, d, dim] >= point[c, d, dim] - barriers[a-1000][b][dim])
            model.addConstr(dif[a, b, c, d, dim] >= - point[c, d, dim] + barriers[a-1000][b][dim])
            model.addConstr(gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] * dist[a, b, c, d])

    L_out = 0
    U_out = 100000

    # p constraints
    for a, b, c, d, e, f in paux.keys():
        if (a, b, c, d) in edges_source:
            L_out = 0
            U_out = 100000

            source = sources[a - 1]
            punto = barriers[c-1000][d]
            L_out = af.estima_L(source, punto)
            U_out = af.estima_U(source, punto)

            # print((L_out, U_out))
            # model.addConstr(paux[a, b, c, d, e, f] <= U_out * aux[a, b, c, d, e, f])
            # model.addConstr(paux[a, b, c, d, e, f] <= dist[a, b, c, d])
            model.addConstr(paux[a, b, c, d, e, f] >= L_out * aux[a, b, c, d, e, f])
            model.addConstr(paux[a, b, c, d, e, f] >= dist[a, b, c, d] - U_out * (1 - aux[a, b, c, d, e, f]))

        if (a, b, c, d) in edges_target:
            L_out = 0
            U_out = 100000

            target = targets[abs(c) - 1]
            punto = barriers[a-1000][b]
            L_out = af.estima_L(target, punto)
            U_out = af.estima_U(target, punto)

            # print((L_out, U_out))
            # model.addConstr(paux[a, b, c, d, e, f] <= U_out * aux[a, b, c, d, e, f])
            # model.addConstr(paux[a, b, c, d, e, f] <= dist[a, b, c, d])
            model.addConstr(paux[a, b, c, d, e, f] >= L_out * aux[a, b, c, d, e, f])
            model.addConstr(paux[a, b, c, d, e, f] >= dist[a, b, c, d] - U_out* (1 -aux[a, b, c, d, e, f]))

        if (a, b, c, d) in edges_barrier:

            L_out = np.linalg.norm(np.array(barriers[a - 1000][b]) - np.array(barriers[c - 1000][d]))
            U_out = np.linalg.norm(np.array(barriers[a - 1000][b]) - np.array(barriers[c - 1000][d]))

            # model.addConstr(paux[a, b, c, d, e, f] <= U_out * aux[a, b, c, d, e, f])
            # model.addConstr(paux[a, b, c, d, e, f] <= dist[a, b, c, d])
            model.addConstr(paux[a, b, c, d, e, f] >= L_out * aux[a, b, c, d, e, f])
            model.addConstr(paux[a, b, c, d, e, f] >= dist[a, b, c, d] - U_out* (1 -aux[a, b, c, d, e, f]))


    model.addConstrs(p[a, b, c, d] == gp.quicksum((abs(e))*paux[a, b, c, d, e, f] for e, f in vertices_target) for a, b, c, d in edges_total)
    # model.addConstrs(p[a, b, c, d] >= x[a, b, c, d]*dist[a, b, c, d] for a, b, c, d in edges_total)
    # flow conservation constraints
    # for index in y_index:
    #     if len(index) == 3:
    # model.addConstr(gp.quicksum(y[tupla] for tupla in E if tupla[0] == -1) == 1)
    #
    # for v, i in VB:
    #     tuplas_salen = gp.quicksum([y[a, b, c, d] for a, b, c, d in E if (a == v and b == i)])
    #     tuplas_entran = gp.quicksum([y[a, b, c, d] for a, b, c, d in E if (c == v and d == i)])
    #
    #     model.addConstr(tuplas_salen - tuplas_entran == 0)
    #
    # model.addConstr(- gp.quicksum(y[a, b, c, d] for a, b, c, d in E if c == -2) == -1)

    # We only open k facilities
    model.addConstr(y.sum('*') == k)

    # All the targets must be associated to one facility
    model.addConstrs(z.sum('*', '*', c, d) == 1 for c, d in vertices_target)

    # If a facility is not open in S, we can not launch capacity from this source.
    model.addConstrs(z[a, b, c, d] <= y[a, b] for a, b, c, d in z_index)

    # Binarying an integer variable
    model.addConstrs(x[a, b, c, d] == gp.quicksum((abs(e))*aux[a, b, c, d, e, f] for e, f in vertices_target) for a, b, c, d in edges_total)
    # model.addConstrs(aux.sum(a, b, c, d, '*', '*') == 1 for a, b, c, d in edges_total)

    for a, b in vertices_total:
        if (a, b) in vertices_source:
            model.addConstr(gp.quicksum(x[a, b, c, d] for c, d in vertices_total if (a, b, c, d) in edges_total) == z.sum(a, b, '*', '*'))
        elif (a, b) in vertices_barrier:
            model.addConstr(gp.quicksum(x[a, b, c, d] for c, d in vertices_total if (a, b, c, d) in edges_total) - gp.quicksum(x[c, d, a, b] for c, d in vertices_total if ((c, d, a, b) in edges_total)) == 0)
        else:
            model.addConstr(gp.quicksum(x[c, d, a, b] for c, d in vertices_total if (c, d, a, b) in edges_total) == 1)

    model.update()

    objective = gp.quicksum(p[index] for index in edges_source + edges_target) + gp.quicksum(0.5*wL*x[index] for index in x.keys())

    for a, b, c, d in edges_barrier:
        # for e, f in vertices_source:
        #     for g, h in vertices_target:
        objective += dist[a, b, c, d]*x[a, b, c, d]

    model.setObjective(objective, GRB.MINIMIZE)

    second_time = time.time()

    time_elapsed = second_time - first_time

    model.update()

    model.Params.Threads = 6
    model.Params.timeLimit = time_limit # - time_elapsed
    # model.Params.LazyConstraints = 1
    model.Params.NumericFocus = 1
    # model.Params.NonConvex = 2

    # model.write('prueba.lp')
    # model.write('prueba.mps')

    model.optimize()

    results = [len(barriers), len(sources), len(barriers), k, wL, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    if init:
        try:
            results[-2] = time_h
            results[-1] = objval_h
        except:
            print('No solution obtained by the heuristic')

    if model.Status == 3:
        model.computeIIS()
        model.write('infeasible_constraints.ilp')
        return results

    if model.SolCount == 0:
        return results

    # model.write('solution.sol')

    results[-6] = model.getAttr('MIPGap')
    results[-5] = model.Runtime + time_elapsed
    results[-4] = model.getAttr('NodeCount')
    results[-3] = model.ObjVal

    x_indices = [(index, x[index].X) for index in x.keys() if x[index].X > 0.5]
    dist_indices = [(index, dist[index].X) for index in x.keys() if x[index].X > 0.5]
    p_indices = [(index, p[index].X) for index in p.keys() if x[index].X > 0.5]
    paux_indices = [(index, paux[index].X) for index in paux.keys() if x[index[0:4]].X > 0.5]
    aux_indices = [(index, aux[index].X) for index in aux.keys() if aux[index].X > 0.5]

    # print("x_indices:")
    # print(x_indices)
    #
    # print("aux_indices:")
    # print(aux_indices)
    #
    # print("paux_indices:")
    # print(paux_indices)
    #
    # print("p_indices:")
    # print(p_indices)
    #
    # print("dist_indices:")
    # print(dist_indices)

    x_indices = [index for index in x.keys() if x[index].X > 0.5]

    y_indices = [index for index in y.keys() if y[index].X > 0.5]

    z_indices = [index for index in z.keys() if z[index].X > 0.5]

    if log:
        print(x_indices)
        print(y_indices)

    # g_indices = []
    #
    # for index in g_index:
    #     if g[index].X > 0.5:
    #         g_indices.append(g[index])
    #
    # if log:
    #     print(g_indices)

    if picture:
        fig, ax = plt.subplots()

        for b in barriers:
            ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c='red')

        for n in sources + targets:
            ax.add_artist(n.artist)

        points_vals = model.getAttr('x', point)
        print(points_vals)

        if log:
            print(points_vals)

        for a, b in y_indices:
            ax.scatter(point[a, b, 0].X, point[a, b, 1].X, s = 30, c = 'black')

        # for c, d in vertices_source:
        #     ax.scatter(point[c, d, 0].X, point[c, d, 1].X, s = 10, c = 'black')

        segments = []

        for a, b, c, d in x_indices:
            if (a, b, c, d) in edges_source:
                segments.append([point[a, b, 0].X, barriers[c-1000][d][0], point[a, b, 1].X, barriers[c-1000][d][1]])

            if (a, b, c, d) in edges_barrier:
                segments.append([barriers[a-1000][b][0], barriers[c-1000][d][0], barriers[a-1000][b][1], barriers[c-1000][d][1]])

            if (a, b, c, d) in edges_target:
                segments.append([barriers[a-1000][b][0], point[c, d, 0].X, barriers[a-1000][b][1], point[c, d, 1].X])


        # print(segments)

        for segment in segments:
            ax.arrow(segment[0], segment[2], segment[1] - segment[0], segment[3] - segment[2], width=0.1,
                     head_width=1, length_includes_head=True, color='black')

        # plt.axis([-5, 105, -5, 105])
        plt.axis([0, 100, 0, 100])

        ax.set_aspect('equal')
        plt.show()

    return results
