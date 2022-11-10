"""
Matheuristic procedure to handle larger instances. Fixing the centers of the neighbourhoods.

"""
import gurobipy as gp
import itertools
import time
from gurobipy import GRB
import numpy as np
from matplotlib.patches import Circle

import estimacion_M as eM
import auxiliar_functions as af
import neighborhood as neigh
from data import *


def matheuristic(barriers, sources_auxiliar, targets_auxiliar, k, wL=50, A4 = True, log=False, picture=False,
                time_limit=7200, init=False):

    print('Iniciando matheuristico')

    S = len(sources_auxiliar)
    T = len(targets_auxiliar)

    sources = af.preprocess_neighborhoods(sources_auxiliar)
    targets = af.preprocess_neighborhoods(targets_auxiliar)

    first_time = time.time()

    # Indices of the vertices identified with the sources
    vertices_source = list(itertools.product(range(1, len(sources) + 1), range(1)))

    # Indices of the vertices identified with the vertices of the barriers
    vertices_barrier = list(itertools.product(range(1000, 1000 + len(barriers)), range(2)))

    # Indices of the vertices identified with the sources
    vertices_target = list(itertools.product(range(-len(targets), 0), range(1)))
    vertices_target = vertices_target[::-1]

    point = {}
    for a, b in vertices_source + vertices_target + vertices_barrier:
        if (a, b) in vertices_source:
            source = sources[a-1]
            if type(source) is neigh.Circle:
                for dim in range(2):
                    point[a, b, dim] = source.center[dim]

            if type(source) is neigh.Poligonal:
                for dim in range(2):
                    point[a, b, dim] = source.V[0][dim]

        if (a, b) in vertices_target:
            target = targets[abs(a)-1]
            if type(target) is neigh.Circle:
                for dim in range(2):
                    point[a, b, dim] = target.center[dim]

            if type(target) is neigh.Poligonal:
                for dim in range(2):
                    point[a, b, dim] = target.V[0][dim]

        if (a, b) in vertices_barrier:
            for dim in range(2):
                point[a, b, dim] = barriers[a-1000][b][dim]

    # Indices of the edges joining a source and a barrier
    edges_source = []

    for (a, b) in vertices_source:
        for c, d in vertices_barrier:
            # Point of the barrier to check if is visible by the neighborhood
            point2 = barriers[c-1000][d]

            # Neighborhood to check if it is visible by the point
            point1 = [point[a, b, 0], point[a, b, 1]]

            barrier = [point1, point2]

            inter = False
            for barrieri in barriers:
                if af.intersect(barrieri, barrier):
                    inter = True
                    break

            if not (inter):
                edges_source.append((a, b, c, d))

    edges_barrier = []
    for v, i in vertices_barrier:
        for w, j in vertices_barrier:
            if v != w:
                barrier = [barriers[v-1000][i], barriers[w-1000][j]]

                inter = False
                for barrieri in barriers:
                    if af.intersect(barrieri, barrier):
                        inter = True
                        break

                if not (inter):
                    edges_barrier.append((v, i, w, j))

    indices_barriers = [(v, 0, v, 1) for v in range(1000, 1000 + len(barriers))]

    edges_target = []

    for (a, b) in vertices_target:
        for c, d in vertices_barrier:
            # Point of the barrier to check if is visible by the neighborhood
            point2 = barriers[c-1000][d]

            # Neighborhood to check if it is visible by the point
            point1 = [point[a, b, 0], point[a, b, 1]]

            barrier = [point1, point2]

            inter = False
            for barrieri in barriers:
                if af.intersect(barrieri, barrier):
                    inter = True
                    break

            if not (inter):
                edges_target.append((c, d, a, b))

    # Including edges joining source neighbourhood and target neighbourhood

    edges_source_target = []

    if not (A4):

        for (a, b) in vertices_source:
            for (c, d) in vertices_target:
                # Point of the barrier to check if is visible by the neighborhood
                point2 = [point[c, d, 0], point[c, d, 1]]

                # Neighborhood to check if it is visible by the point
                point1 = [point[a, b, 0], point[a, b, 1]]

                barrier = [point1, point2]

                inter = False
                for barrieri in barriers:
                    if af.intersect(barrieri, barrier):
                        inter = True
                        break

                if not (inter):
                    edges_source_target.append((a, b, c, d))

    vertices_total = vertices_source + vertices_barrier + vertices_target
    edges_total = edges_source + edges_barrier + edges_target + edges_source_target

    if log:
        print("vertices_source = " + str(vertices_source))
        print("vertices_barrier = " + str(vertices_barrier))
        print("vertices_target = " + str(vertices_target))

        print("edges_source = " + str(edges_source))
        print("edges_barrier = " + str(edges_barrier))
        print("edges_free = " + str(edges_target))
        print("edges_total = " + str(edges_total))

    epsilon_index = []  # epsilon(S / T, B, i) = 1 si (P_{S/T}, P_B^i)\in E_{S/T}

    # print(epsilon_index)

    # point_index = []
    # for a, b in vertices_source+vertices_target:
    #     for dim in range(2):
    #         point_index.append((a, b, dim))

    x_index = edges_total

    aux_index = []

    for a, b, c, d in edges_total:
        for e, f in vertices_target:
            aux_index.append((a, b, c, d, e, f))


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

    # for index in vertices_neighborhood:
    #     for dim in range(2):
    #         p_index.append((index[0], index[1], dim))



    dist = {}

    for a, b, c, d in edges_total:

        dist[a, b, c, d] = np.linalg.norm(np.array([point[a, b, 0], point[a, b, 1]]) - np.array([point[c, d, 0], point[c, d, 1]]))

        # if (a, b, c, d) in edges_source:
        #     source = sources[a-1]
        #
        #     if type(source) is neigh.Circle:
        #         dist[a, b, c, d] = np.linalg.norm(np.array(source.center)-np.array(barriers[c-1000][d]))
        #     if type(source) is neigh.Poligonal:
        #         dist[a, b, c, d] = np.linalg.norm(np.array(source.V[0])-np.array(barriers[c-1000][d]))
        #
        # elif (a, b, c, d) in edges_barrier:
        #     dist[a, b, c, d] = np.linalg.norm(np.array(barriers[a-1000][b]) - np.array(barriers[c - 1000][d]))
        #
        # elif (a, b, c, d) in edges_target:
        #     target = targets[abs(c)-1]
        #
        #     if type(target) is neigh.Circle:
        #         dist[a, b, c, d] = np.linalg.norm(np.array(target.center)-np.array(barriers[a-1000][b]))
        #     if type(target) is neigh.Poligonal:
        #         dist[a, b, c, d] = np.linalg.norm(np.array(target.V[0])-np.array(barriers[a-1000][b]))
        #
        # elif (a, b, c, d) in edges_source_target:
        #     source = sources[a-1]
        #     target = targets[abs(c) - 1]
        #
        #     if type(source) is neigh.Circle and type(target) is neigh.Circle:
        #         dist[a, b, c, d] = np.linalg.norm(np.array(source.center) - np.array(target.center))
        #     elif type(source) is neigh.Circle and type(target) is neigh.Poligonal:
        #         dist[a, b, c, d] = np.linalg.norm(np.array(source.center) - np.array(target.V[0]))
        #     elif type(source) is neigh.Poligonal and type(target) is neigh.Circle:
        #         dist[a, b, c, d] = np.linalg.norm(np.array(source.V[0]) - np.array(target.center))
        #     elif type(source) is neigh.Poligonal and type(target) is neigh.Poligonal:
        #         dist[a, b, c, d] = np.linalg.norm(np.array(source.V[0]) - np.array(target.V[0]))


    model = gp.Model('Model: H-KMedian-N')

    x = model.addVars(x_index, vtype=GRB.INTEGER, lb=0.0, ub=len(targets), name='x')
    aux = model.addVars(aux_index, vtype=GRB.BINARY, name='aux')
    y = model.addVars(y_index, vtype=GRB.BINARY, name='y')
    z = model.addVars(z_index, vtype=GRB.BINARY, name='z')

    model.update()

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
            model.addConstr(x.sum(a, b, '*', '*') == z.sum(a, b, '*', '*'))
        elif (a, b) in vertices_barrier:
            model.addConstr(x.sum(a, b, '*', '*') - x.sum('*', '*', a, b) == 0)
        else:
            model.addConstr(x.sum('*', '*', a, b) == 1)

    model.update()

    objective = gp.quicksum(dist[index]*x[index] for index in edges_total) + gp.quicksum(0.5*wL*x[index] for index in x.keys())

    # for a, b, c, d in edges_barrier:
    #     # for e, f in vertices_source:
    #     #     for g, h in vertices_target:
    #     objective += dist[a, b, c, d]*x[a, b, c, d]

    model.setObjective(objective, GRB.MINIMIZE)

    model.update()

    model.Params.Threads = 6
    model.Params.timeLimit = time_limit  # - time_elapsed
    model.Params.NumericFocus = 1

    if not (A4):
        model.Params.NonConvex = 2


    model.optimize()

    second_time = time.time()

    time_elapsed = second_time - first_time

    model.write('solution.sol')
    results = [np.nan, np.nan]

    if model.Status == 3:
        model.computeIIS()
        model.write('infeasible_constraints.ilp')
        return results

    if model.SolCount == 0:
        return results

    # model.write('solution.sol')

    time_h = model.getAttr('Runtime')
    objval_h = model.ObjVal

    x_indices = [index for index in x.keys() if x[index].X > 0.5]

    y_indices = [index for index in y.keys() if y[index].X > 0.5]

    z_indices = [index for index in z.keys() if z[index].X > 0.5]

    if log:
        print(x_indices)
        print(y_indices)
        print(z_indices)

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

        for a, b in y_indices:
            ax.scatter(point[a, b, 0], point[a, b, 1], s = 30, c = 'black')

        # for c, d in vertices_source:
        #     ax.scatter(point[c, d, 0].X, point[c, d, 1].X, s = 10, c = 'black')

        segments = []

        for a, b, c, d in x_indices:
            if (a, b, c, d) in edges_source:
                segments.append([point[a, b, 0], barriers[c-1000][d][0], point[a, b, 1], barriers[c-1000][d][1]])

            if (a, b, c, d) in edges_barrier:
                segments.append([barriers[a-1000][b][0], barriers[c-1000][d][0], barriers[a-1000][b][1], barriers[c-1000][d][1]])

            if (a, b, c, d) in edges_target:
                segments.append([barriers[a-1000][b][0], point[c, d, 0], barriers[a-1000][b][1], point[c, d, 1]])


        # print(segments)

        for segment in segments:
            ax.arrow(segment[0], segment[2], segment[1] - segment[0], segment[3] - segment[2], width=0.1,
                     head_width=1, length_includes_head=True, color='black')

        # plt.axis([-5, 105, -5, 105])
        plt.axis([0, 100, 0, 100])

        ax.set_aspect('equal')
        plt.show()

    return time_h, objval_h, x_indices, y_indices, z_indices


