# TSPN-B

import gurobipy as gp
import itertools
import time
from gurobipy import GRB
from matplotlib.patches import Circle

import estimacion_M as eM
import auxiliar_functions as af
import neighborhood as neigh
from data import *
from matheuristic import matheuristic


def h_kmedian_n(barriers, sources, targets, k, single = False, wL=50, lazy=True, A4=True, prepro=True, log=False, picture=False, time_limit=7200, init=False):
    
    # Size of the set of sources and targets
    S = len(sources)
    T = len(targets)

    first_time = time.time()

    # Indices of the vertices identified with the sources
    vertices_source = list(itertools.product(range(1, len(sources) + 1), range(1)))

    # Indices of the vertices identified with the vertices of the barriers
    vertices_barrier = list(itertools.product(range(1000, 1000 + len(barriers)), range(2)))

    # Indices of the vertices identified with the targets
    vertices_target = list(itertools.product(range(-len(targets), 0), range(1)))
    vertices_target = vertices_target[::-1]

    # print('Number of sources: ' + str(S))
    # print('Number of targets: ' + str(S))
    # print('Number of barriers: ' + str(len(barriers)))

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
            else:
                edges_source.append((a, b, c, d))
    
    second_time = time.time() - first_time
    print('\nTime_source_barrier: ' + str(second_time))

    start_barrier_barrier = time.time()
    edges_barrier = []
    for v, i in vertices_barrier:
        for w, j in vertices_barrier:
            if v > w:
                barrier = [barriers[v-1000][i], barriers[w-1000][j]]

                if np.linalg.norm(np.array(barriers[v-1000][i]) - np.array(barriers[w-1000][j])) >= 0.5:

                    intersect = False
                    for barrieri in barriers:
                        if af.intersect(barrieri, barrier):
                            intersect = True
                            break

                    if not (intersect):
                        edges_barrier.append((v, i, w, j))
                        edges_barrier.append((w, j, v, i))
                    else:
                        pass
    
    third_time = time.time() - start_barrier_barrier
    print('Time_barrier_barrier: ' + str(third_time))

    indices_barriers = [(v, 0, v, 1) for v in range(1000, 1000 + len(barriers))]

    edges_target = []

    start_target_barrier = time.time()
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

    # print(edges_target)

    fourth_time = time.time() - start_target_barrier
    print('Time_barrier_target: ' + str(fourth_time))

    start_source_target = time.time()
    # Including edges joining source neighbourhood and target neighbourhood
    edges_source_target = []

    if not(A4):
        for (a, b) in vertices_source:
            for (c, d) in vertices_target:
                neighborhood1 = sources[a - 1]
                neighborhood2 = targets[abs(c) - 1]

                if af.canseeN(neighborhood1, neighborhood2, barriers):
                    edges_source_target.append((a, b, c, d))

    fifth_time = time.time() - start_source_target
    print('Time_source_target: ' + str(fifth_time))

    vertices_total = vertices_source + vertices_barrier + vertices_target

    if not(A4):
        edges_total = edges_source + edges_barrier + edges_target + edges_source_target
    else:
        edges_total = edges_source + edges_barrier + edges_target

    if log:
        print("vertices_source = " + str(vertices_source))
        print("vertices_barrier = " + str(vertices_barrier))
        print("vertices_target = " + str(vertices_target))

        print("edges_source = " + str(edges_source))
        print("edges_barrier = " + str(edges_barrier))
        print("edges_free = " + str(edges_target))
        print("edges_total = " + str(edges_total))

    # Indices of the points
    point_index = []
    for a, b in vertices_source + vertices_target:
        for dim in range(2):
            point_index.append((a, b, dim))

    shape_barriers = np.array(barriers).shape[0]
    bounds_min, bounds_max = min([0, np.array(barriers).reshape(shape_barriers*2, 2).min()]), max([100, np.array(barriers).reshape(shape_barriers*2, 2).max()])


    # Indices of alpha variables
    alpha_index = []

    for a, b in vertices_barrier:
        for c, d, e, f in edges_source + edges_target + edges_source_target:
            alpha_index.append((a, b, c, d, e, f))
    for a, b in vertices_source + vertices_target:
        for (c, d, e, f) in indices_barriers:
            alpha_index.append((a, b, c, d, e, f))

    print("\nCardinality of alpha set: " + str(len(alpha_index)))

    if log:
        print("alpha = " + str(alpha_index))

    # Indices of beta variables
    beta_index = []

    for a, b, c, d in edges_source + edges_target + edges_source_target:
        for e, f, g, h in indices_barriers:
            beta_index.append((a, b, c, d, e, f, g, h))
            beta_index.append((e, f, g, h, a, b, c, d))

    print("Cardinality of beta set: " + str(len(beta_index)))

    # Indices of gamma variables
    gamma_index = beta_index

    if log:
        print("beta/gamma = " + str(beta_index))

    # Indices of delta variables
    delta_index = []

    for a, b, c, d in edges_source + edges_target + edges_source_target:
        for e, f, g, h, in indices_barriers:
            delta_index.append((a, b, c, d, e, f, g, h))

    print("Cardinality of delta set: " + str(len(delta_index)))

    if log:
        print("delta = " + str(delta_index))

    # Indices of epsilon variables
    epsilon_index = []

    for a, b, c, d in edges_source + edges_target + edges_source_target:
        epsilon_index.append((a, b, c, d))

    print("Cardinality of epsilon set: " + str(len(epsilon_index)))

    # Assignment variables
    x_index = []

    for a, b in vertices_source:
        for c, d in vertices_target:
            x_index.append((a, b, c, d))

    # Allocating variables
    y_index = vertices_source

    # Flow variables
    flow_index = []
    aux_index = []

    # Single commodity flow case
    if single:
        for a, b, c, d in edges_total:
            flow_index.append((a, b, c, d))
            for e, f in vertices_target:
                aux_index.append((a, b, c, d, e, f))
        
        p_index = edges_total

    else:
    # Multi commodity flow case
        # for a, b, c, d in edges_total:
        #     for e, f, g, h in x_index:
        #         if ((a, b, c, d) in edges_source and a == e) or ((a, b, c, d) in edges_target and c == g) or ((a, b, c, d) in edges_barrier) or ((a, b, c, d) in edges_source_target and a == e and c == g):
        #             flow_index.append((a, b, c, d, e, f, g, h))

        for a, b, c, d in edges_source:
            for g, h in vertices_target:
                flow_index.append((a, b, c, d, a, 0, g, h))
        
        for a, b, c, d in edges_target:
            for e, f in vertices_source:
                flow_index.append((a, b, c, d, e, f, c, 0))
        
        for a, b, c, d in edges_barrier:
            for e, f, g, h in x_index:
                flow_index.append((a, b, c, d, e, f, g, h))
        
        for a, b, c, d in edges_source_target:
            flow_index.append((a, b, c, d, a, 0, c, 0))

        aux_index = flow_index

        p_index = x_index

    paux_index = aux_index

    if log:
        print("x_index = " + str(x_index))
        print("y_index = " + str(y_index))
        print("flow_index = " + str(flow_index))


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

    if lazy:
        def elimcuts(model, where):
            if where == GRB.Callback.MIPSOL:
                flows = model.cbGetSolution(model._flow)

                alphas = model.cbGetSolution(model._alpha)
                betas = model.cbGetSolution(model._beta)
                gammas = model.cbGetSolution(model._gamma)

                points = model.cbGetSolution(model._point)

                indices = []

                # Alpha constraint
                L = -100000
                U = 100000

                # Loop to add constraints
                flow_indices = [index for index in flows.keys() if flows[index] > 0.5]

                for a, b, c, d in flow_indices:
                    if (a, b, c, d) in edges_source:
                        segment = [[points[a, b, 0], points[a, b, 1]], barriers[c - 1000][d]]



                        for e in range(1000, 1000 + len(barriers)):
                            det1 = af.determinant([points[a, b, 0], points[a, b, 1]], barriers[e - 1000][0], barriers[e-1000][1])
                            det2 = af.determinant(barriers[c-1000][d], barriers[e - 1000][0], barriers[e-1000][1])

                            det3 = af.determinant(barriers[e-1000][0], [points[a, b, 0], points[a, b, 1]], barriers[c-1000][d])
                            det4 = af.determinant(barriers[e-1000][1], [points[a, b, 0], points[a, b, 1]], barriers[c-1000][d])

                            if det1*det2 < 0 and det3*det4 < 0:
                                # print(af.determinant([model._point[a, b, 0], model._point[a, b, 1]], barriers[e - 1000][0], barriers[e - 1000][1]))
                                L1, U1 = eM.estima_M_alpha1(sources[a - 1], barriers[e - 1000][0], barriers[e - 1000][1])
                                model.cbLazy(
                                    (1 - model._alpha[a, b, e, 0, e, 1]) * L1 <= af.determinant([model._point[a, b, 0], model._point[a, b, 1]], barriers[e - 1000][0], barriers[e - 1000][1]))
                                model.cbLazy(
                                    -U1 * model._alpha[a, b, e, 0, e, 1] <= -af.determinant([model._point[a, b, 0], model._point[a, b, 1]], barriers[e - 1000][0], barriers[e - 1000][1]))

                                for f in range(2):
                                    L, U = eM.estima_M_alpha2(barriers[e - 1000][f], sources[a - 1], barriers[c - 1000][d])

                                    model.cbLazy(
                                        (1 - model._alpha[e, f, a, b, c, d]) * L <= af.determinant(barriers[e - 1000][f], [model._point[a, b, 0], model._point[a, b, 1]], barriers[c - 1000][d]))
                                    model.cbLazy(
                                        - U * model._alpha[e, f, a, b, c, d] <= - af.determinant(barriers[e - 1000][f], [model._point[a, b, 0], model._point[a, b, 1]], barriers[c - 1000][d]))


                    elif (a, b, c, d) in edges_target:
                        segment = [barriers[a - 1000][b], [points[c, d, 0], points[c, d, 1]]]

                        for e in range(1000, 1000 + len(barriers)):
                            det1 = af.determinant(barriers[a - 1000][b], barriers[e - 1000][0], barriers[e - 1000][1])
                            det2 = af.determinant([points[c, d, 0], points[c, d, 1]], barriers[e - 1000][0], barriers[e - 1000][1])

                            det3 = af.determinant(barriers[e - 1000][0], barriers[a - 1000][b], [points[c, d, 0], points[c, d, 1]])
                            det4 = af.determinant(barriers[e - 1000][1], barriers[a - 1000][b], [points[c, d, 0], points[c, d, 1]])

                            if det1 * det2 < 0 and det3 * det4 < 0:
                                L1, U1 = eM.estima_M_alpha1(targets[abs(c) - 1], barriers[e - 1000][0], barriers[e - 1000][1])

                                model.cbLazy(
                                    (1 - model._alpha[c, d, e, 0, e, 1]) * L1 <= af.determinant([model._point[c, d, 0], model._point[c, d, 1]], barriers[e - 1000][0], barriers[e - 1000][1]))
                                model.cbLazy(
                                    -U1 * model._alpha[c, d, e, 0, e, 1] <= -af.determinant([model._point[c, d, 0], model._point[c, d, 1]], barriers[e - 1000][0], barriers[e - 1000][1]))

                                for f in range(2):
                                    L, U = eM.estima_M_alpha3(barriers[e - 1000][f], barriers[a - 1000][b], targets[abs(c) - 1])
                                    model.cbLazy(
                                        (1 - model._alpha[e, f, a, b, c, d]) * L <= af.determinant(barriers[e - 1000][f], barriers[a - 1000][b], [model._point[c, d, 0], model._point[c, d, 1]]))
                                    model.cbLazy(
                                        - U * model._alpha[e, f, a, b, c, d] <= - af.determinant(barriers[e - 1000][f], barriers[a - 1000][b], [model._point[c, d, 0], model._point[c, d, 1]]))

                            # print("Llega hasta aqui")

                    elif (a, b, c, d) in edges_source_target:
                        segment = [[points[a, b, 0], points[a, b, 1]], [points[c, d, 0], points[c, d, 1]]]

                        for barrier in barriers:
                            if af.intersect(segment, barrier):
                                indices.append((a, b, c, d))
                                model.cbLazy(model._x[a, b, c, d] <= 0.5)

                    model.cbLazy(
                        model._beta[a, b, c, d, e, 0, e, 1] == 2 * model._gamma[a, b, c, d, e, 0, e, 1] - model._alpha[
                            a, b, e, 0, e, 1] - model._alpha[c, d, e, 0, e, 1] + 1)
                    model.cbLazy(
                        model._beta[e, 0, e, 1, a, b, c, d] == 2 * model._gamma[e, 0, e, 1, a, b, c, d] - model._alpha[
                            e, 0, a, b, c, d] - model._alpha[e, 1, a, b, c, d] + 1)

                    model.cbLazy(model._gamma[a, b, c, d, e, 0, e, 1] <= model._alpha[a, b, e, 0, e, 1])
                    model.cbLazy(model._gamma[a, b, c, d, e, 0, e, 1] <= model._alpha[c, d, e, 0, e, 1])
                    model.cbLazy(model._gamma[a, b, c, d, e, 0, e, 1] >= model._alpha[a, b, e, 0, e, 1] +
                                 model._alpha[c, d, e, 0, e, 1] - 1)

                    model.cbLazy(model._gamma[e, 0, e, 1, a, b, c, d] <= model._alpha[e, 0, a, b, c, d])
                    model.cbLazy(model._gamma[e, 0, e, 1, a, b, c, d] <= model._alpha[e, 1, a, b, c, d])
                    model.cbLazy(model._gamma[e, 0, e, 1, a, b, c, d] >= model._alpha[e, 0, a, b, c, d] +
                                 model._alpha[e, 1, a, b, c, d] - 1)

                    model.cbLazy(
                        model._beta[a, b, c, d, e, 0, e, 1] + model._beta[e, 0, e, 1, a, b, c, d] >= 1)


    def callback(model, where):
        if where == GRB.Callback.MIPSOL:
            if model.cbGet(GRB.Callback.MIPSOL_SOLCNT) == 0:
                # creates new model attribute '_startobjval'
                model._startobjval = model.cbGet(GRB.Callback.MIPSOL_OBJ)
                model._starttime = model.cbGet(GRB.Callback.RUNTIME)

                # model.terminate()

    model = gp.Model('Model: H-KMedian-N')

    p = model.addVars(p_index, vtype=GRB.CONTINUOUS, lb=0.0, name='p')
    paux = model.addVars(paux_index, vtype=GRB.CONTINUOUS, lb=0.0, name='paux')
    x = model.addVars(x_index, vtype=GRB.BINARY, name='x')
    aux = model.addVars(aux_index, vtype=GRB.BINARY, name='aux')

    y = model.addVars(y_index, vtype=GRB.BINARY, name='y')
    if single:
        flow = model.addVars(flow_index, vtype=GRB.INTEGER, lb=0.0, ub=T, name='flow')
    else:
        flow = model.addVars(flow_index, vtype=GRB.BINARY, name='flow')
    dist = model.addVars(dist_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dist')
    dif = model.addVars(dif_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif')

    alpha = model.addVars(alpha_index, vtype=GRB.BINARY, name='alpha')
    beta = model.addVars(beta_index, vtype=GRB.BINARY, name='beta')
    gamma = model.addVars(gamma_index, vtype=GRB.BINARY, name='gamma')

    if not (lazy):
        delta = model.addVars(delta_index, vtype=GRB.BINARY, name='delta')
        epsilon = model.addVars(epsilon_index, vtype=GRB.BINARY, name='epsilon')

    point = model.addVars(point_index, vtype=GRB.CONTINUOUS, lb=min([0, bounds_min]), name='point')

    d_inside = model.addVars(d_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='d_inside')
    dif_inside = model.addVars(dif_inside_index, vtype=GRB.CONTINUOUS, lb=0.0, name='dif_inside')

    model.update()

    if init:
        # Solving the problem by taking the center of the neighbourhoods.
        time_h, objval_h, x_start, y_start, flow_start, point_vals = matheuristic(barriers, sources, targets, edges_source, edges_target, edges_barrier, edges_source_target, k, single = single, wL = wL, A4=A4, time_limit = 100, picture = False)

        # model.read('solution.sol')
        # print(flow_start)

        # print(x_start)
        for index in x_start:
            x[index].start = 1
        
        for index in y_start:
            y[index].start = 1
        
        for index in flow_start:
            flow[index].start = 1
        
        for index in vertices_source + vertices_target:
            for dim in range(2):
                point[index[0], index[1], dim].start = point_vals[index[0], index[1], dim]

    L = -10000
    U = 10000

    if not (lazy):
        # alpha constraints
        print('\nSetting alpha constraints: ')

        start_alpha = time.time()

        alpha_counter = 0
        
        for a, b in vertices_source:
            for (c, d, e, f) in indices_barriers:
                L, U = eM.estima_M_alpha1(sources[a-1], barriers[c-1000][d], barriers[e-1000][f])

                model.addConstr(
                    af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c - 1000][d], barriers[e - 1000][f]) >= (1 - alpha[a, b, c, d, e, f]) * L)

                model.addConstr(
                    af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c - 1000][d], barriers[e - 1000][f]) <= U * alpha[a, b, c, d, e, f])

                alpha_counter += 2

        for a, b in vertices_target:
            for (c, d, e, f) in indices_barriers:
                L, U = eM.estima_M_alpha1(targets[abs(a)-1], barriers[c-1000][d], barriers[e-1000][f])

                model.addConstr(
                    af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c - 1000][d], barriers[e - 1000][f]) >= (1 - alpha[a, b, c, d, e, f]) * L)
                
                model.addConstr(
                    af.determinant([point[a, b, 0], point[a, b, 1]], barriers[c - 1000][d], barriers[e - 1000][f]) <= U * alpha[a, b, c, d, e, f])

                alpha_counter += 2

        for a, b in vertices_barrier:
            for c, d, e, f in edges_source:
                L, U = eM.estima_M_alpha2(barriers[a-1000][b], sources[c-1], barriers[e-1000][f])
                
                model.addConstr(
                    af.determinant(barriers[a - 1000][b], [point[c, d, 0], point[c, d, 1]], barriers[e - 1000][f]) >= (1 - alpha[a, b, c, d, e, f]) * L)
                
                model.addConstr(
                    af.determinant(barriers[a - 1000][b], [point[c, d, 0], point[c, d, 1]], barriers[e - 1000][f]) <= U * alpha[a, b, c, d, e, f])
                
                alpha_counter += 2

        for a, b in vertices_barrier:
            for c, d, e, f in edges_target:
                L, U = eM.estima_M_alpha3(barriers[a-1000][b], barriers[c-1000][d], targets[abs(e)-1])
                
                model.addConstr(
                    af.determinant(barriers[a - 1000][b], barriers[c - 1000][d], [point[e, f, 0], point[e, f, 1]]) >= (1 - alpha[a, b, c, d, e, f]) * L)
                
                model.addConstr(
                    af.determinant(barriers[a - 1000][b], barriers[c - 1000][d], [point[e, f, 0], point[e, f, 1]]) <= alpha[a, b, c, d, e, f] * U)

                alpha_counter += 2

        for a, b in vertices_barrier:
            for c, d, e, f in edges_source_target:                
                L, U = eM.estima_M_alpha4(barriers[a-1000][b], sources[c-1], targets[abs(e)-1])
                
                model.addConstr(
                    af.determinant(barriers[a - 1000][b], [point[c, d, 0], point[c, d, 1]], [point[e, f, 0], point[e, f, 1]]) >= (1 - alpha[a, b, c, d, e, f]) * L)
                
                model.addConstr(
                    af.determinant(barriers[a - 1000][b], [point[c, d, 0], point[c, d, 1]], [point[e, f, 0], point[e, f, 1]]) <=  alpha[a, b, c, d, e, f] * U)

                alpha_counter += 2

        print('Spent time: ' + str(time.time() - start_alpha))
        print('Number of alpha constraints: ' + str(alpha_counter))

        model.update()

        print('\nSetting beta and gamma constraints: ')

        start_beta = time.time()
        
        beta_counter = 0
        gamma_counter = 0

        # beta constraints
        for a, b, c, d in edges_source:
            for e, f, g, h in indices_barriers:
                value = (af.determinant(barriers[c - 1000][d], barriers[e - 1000][f], barriers[g - 1000][h]) > 0)*1

                model.addConstr(
                    beta[a, b, c, d, e, f, g, h] == 2 * alpha[a, b, e, f, g, h] * value - alpha[a, b, e, f, g, h] - value + 1)

                model.addConstr(
                    beta[e, f, g, h, a, b, c, d] == 2 * gamma[e, f, g, h, a, b, c, d] - alpha[e, f, a, b, c, d] - alpha[g, h, a, b, c, d] + 1)

                beta_counter += 2

                model.addConstr(gamma[e, f, g, h, a, b, c, d] <= alpha[e, f, a, b, c, d])
                model.addConstr(gamma[e, f, g, h, a, b, c, d] <= alpha[g, h, a, b, c, d])
                model.addConstr(gamma[e, f, g, h, a, b, c, d] >= alpha[e, f, a, b, c, d] + alpha[g, h, a, b, c, d] - 1)

                gamma_counter += 3

        for a, b, c, d in edges_target:
            for e, f, g, h in indices_barriers:
                value = (af.determinant(barriers[a - 1000][b], barriers[e - 1000][f], barriers[g - 1000][h]) > 0)*1

                model.addConstr(
                    beta[a, b, c, d, e, f, g, h] == 2 * value * alpha[c, d, e, f, g, h] - value - alpha[c, d, e, f, g, h] + 1)

                model.addConstr(
                    beta[e, f, g, h, a, b, c, d] == 2 * gamma[e, f, g, h, a, b, c, d] - alpha[e, f, a, b, c, d] - alpha[g, h, a, b, c, d] + 1)

                beta_counter += 2

                model.addConstr(gamma[e, f, g, h, a, b, c, d] <= alpha[e, f, a, b, c, d])
                model.addConstr(gamma[e, f, g, h, a, b, c, d] <= alpha[g, h, a, b, c, d])
                model.addConstr(gamma[e, f, g, h, a, b, c, d] >= alpha[e, f, a, b, c, d] + alpha[g, h, a, b, c, d] - 1)

                gamma_counter += 3

        for a, b, c, d in edges_source_target:
            for e, f, g, h in indices_barriers:
                model.addConstr(
                    beta[a, b, c, d, e, f, g, h] == 2 * gamma[a, b, c, d, e, f, g, h] - alpha[a, b, e, f, g, h]  - alpha[c, d, e, f, g, h] + 1)

                model.addConstr(gamma[a, b, c, d, e, f, g, h] <= alpha[a, b, e, f, g, h])
                model.addConstr(gamma[a, b, c, d, e, f, g, h] <= alpha[c, d, e, f, g, h])
                model.addConstr(gamma[a, b, c, d, e, f, g, h] >= alpha[a, b, e, f, g, h] + alpha[c, d, e, f, g, h] - 1)

                model.addConstr(
                    beta[e, f, g, h, a, b, c, d] == 2 * gamma[e, f, g, h, a, b, c, d] - alpha[e, f, a, b, c, d] - alpha[g, h, a, b, c, d] + 1)

                model.addConstr(gamma[e, f, g, h, a, b, c, d] <= alpha[e, f, a, b, c, d])
                model.addConstr(gamma[e, f, g, h, a, b, c, d] <= alpha[g, h, a, b, c, d])
                model.addConstr(gamma[e, f, g, h, a, b, c, d] >= alpha[e, f, a, b, c, d] + alpha[g, h, a, b, c, d] - 1)

                beta_counter += 2
                gamma_counter += 3

        print('Spent time: ' + str(time.time() - start_beta))
        print('Number of beta constraints: ' + str(beta_counter))
        print('Number of gamma constraints: ' + str(gamma_counter))
        
        print('\nSetting delta constraints: ')
        
        start_delta = time.time()
        delta_counter = 0

        # delta constraints
        for a, b, c, d in edges_source + edges_target + edges_source_target:
            for e, f, g, h in indices_barriers:
                model.addConstr(
                    0.5 * (beta[a, b, c, d, e, f, g, h] + beta[e, f, g, h, a, b, c, d]) <= delta[a, b, c, d, e, f, g, h])
                model.addConstr(
                    2 * (beta[a, b, c, d, e, f, g, h] + beta[e, f, g, h, a, b, c, d]) >= delta[a, b, c, d, e, f, g, h])
                delta_counter += 2

        print('Spent time: ' + str(time.time() - start_delta))
        print('Number of delta constraints: ' + str(delta_counter))

        # epsilon constraints

        print('\nSetting epsilon constraints: ')
        
        start_epsilon = time.time()
        epsilon_counter = 0

        for a, b, c, d in epsilon_index:
            model.addConstr(delta.sum(a, b, c, d, '*', '*', '*', '*') - len(barriers) + 1 <= epsilon[a, b, c, d])
            model.addConstr(len(barriers) * epsilon[a, b, c, d] <= delta.sum(a, b, c, d, '*', '*', '*', '*'))

            epsilon_counter += 2

            if single:
                model.addConstr(flow[a, b, c, d] <= T*epsilon[a, b, c, d])

                epsilon_counter += 1

        if not(single):
            for a, b, c, d in edges_source:
                for g, h in vertices_target:
                    model.addConstr(flow[a, b, c, d, a, 0, g, h] <= epsilon[a, b, c, d])

                    epsilon_counter += 1
                    
            for a, b, c, d in edges_target:
                for e, f in vertices_source:
                    model.addConstr(flow[a, b, c, d, e, f, c, 0] <= epsilon[a, b, c, d])
                    
                    epsilon_counter += 1
            
            for a, b, c, d in edges_source_target:
                model.addConstr(flow[a, b, c, d, a, 0, c, 0] <= epsilon[a, b, c, d])
                
                epsilon_counter += 1

        print('Spent time: ' + str(time.time() - start_epsilon))
        print('Number of epsilon constraints: ' + str(epsilon_counter))       

        model.update()

    print('\nSetting neighbourhood constraints: ')
        
    start_neighbourhood = time.time()

    # neighbourhood constraints
    for a, b, dim in point.keys():
        if (a, b) in vertices_source:
            source = sources[a - 1]

            if type(source) is neigh.Circle:
                model.addConstr(dif_inside[a, b, dim] >= point[a, b, dim] - source.center[dim])
                model.addConstr(dif_inside[a, b, dim] >= source.center[dim] - point[a, b, dim])
                model.addConstr(gp.quicksum(dif_inside[a, b, dim] * dif_inside[a, b, dim] for dim in range(2)) <= d_inside[a, b] * d_inside[a, b])
                model.addConstr(d_inside[a, b] <= source.radii)

            if type(source) is neigh.Poligonal:
                model.addConstrs(point[a, b, dim] == landa[a, b] * source.V[0][dim] + (1 - landa[a, b]) * source.V[1][dim] for dim in range(2))

        elif (a, b) in vertices_target:
            target = targets[abs(a) - 1]

            if type(target) is neigh.Circle:
                model.addConstr(dif_inside[a, b, dim] >= point[a, b, dim] - target.center[dim])
                model.addConstr(dif_inside[a, b, dim] >= target.center[dim] - point[a, b, dim])
                model.addConstr(gp.quicksum(dif_inside[a, b, dim] * dif_inside[a, b, dim] for dim in range(2)) <= d_inside[a, b] * d_inside[a, b])
                model.addConstr(d_inside[a, b] <= target.radii)

            if type(target) is neigh.Poligonal:
                model.addConstrs(point[a, b, dim] == landa[a, b] * target.V[0][dim] + (1 - landa[a, b]) * target.V[1][dim] for dim in range(2))

    print('Spent time: ' + str(time.time() - start_neighbourhood))

    print('\nSetting distance constraints: ')
    
    start_distance = time.time()
    # distance constraints constraints
    for a, b, c, d, dim in dif_index:
        if (a, b, c, d) in edges_barrier:
            dist[a, b, c, d] = np.linalg.norm(np.array(barriers[a-1000][b]) - np.array(barriers[c-1000][d]))

        elif (a, b, c, d) in edges_source:
            model.addConstr(dif[a, b, c, d, dim] >= point[a, b, dim] - barriers[c - 1000][d][dim])
            model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + barriers[c - 1000][d][dim])
            model.addConstr(gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] * dist[a, b, c, d])

        elif (a, b, c, d) in edges_target:
            model.addConstr(dif[a, b, c, d, dim] >= point[c, d, dim] - barriers[a - 1000][b][dim])
            model.addConstr(dif[a, b, c, d, dim] >= - point[c, d, dim] + barriers[a - 1000][b][dim])
            model.addConstr(gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] * dist[a, b, c, d])

        else:
            model.addConstr(dif[a, b, c, d, dim] >= point[a, b, dim] - point[c, d, dim])
            model.addConstr(dif[a, b, c, d, dim] >= - point[a, b, dim] + point[c, d, dim])
            model.addConstr(gp.quicksum(dif[a, b, c, d, dim] * dif[a, b, c, d, dim] for dim in range(2)) <= dist[a, b, c, d] * dist[a, b, c, d])

    print('Spent time: ' + str(time.time() - start_distance))

    L_out = 0
    U_out = 100000

    print('\nSetting product constraints: ')
    start_product = time.time()

    # product constraints
    if single:
        for a, b, c, d in edges_source:
            source = sources[a - 1]
            punto = neigh.Punto(barriers[c-1000][d])
            L_out, U_out = eM.estima_M_complete(source, punto)
            for e, f in vertices_target:
                model.addConstr(paux[a, b, c, d, e, f] >= L_out * aux[a, b, c, d, e, f])
                model.addConstr(paux[a, b, c, d, e, f] >= dist[a, b, c, d] - U_out * (1 - aux[a, b, c, d, e, f]))

        for a, b, c, d in edges_target:
            target = targets[abs(c) - 1]
            punto = neigh.Punto(barriers[a-1000][b])
            L_out, U_out = eM.estima_M_complete(target, punto)
            for e, f in vertices_target:
                model.addConstr(paux[a, b, c, d, e, f] >= L_out * aux[a, b, c, d, e, f])
                model.addConstr(paux[a, b, c, d, e, f] >= dist[a, b, c, d] - U_out* (1 -aux[a, b, c, d, e, f]))

        for a, b, c, d in edges_barrier:
            # L_out, U_out = eM.estima_M_complete(neigh.Punto(np.array(barriers[a - 1000][b])), neigh.Punto(np.array(barriers[c - 1000][d])))

            for e, f in vertices_target:
                # model.addConstr(paux[a, b, c, d, e, f] >= L_out * aux[a, b, c, d, e, f])
                # model.addConstr(paux[a, b, c, d, e, f] >= dist[a, b, c, d] - U_out * (1 - aux[a, b, c, d, e, f]))
                paux[a, b, c, d, e, f] = dist[a, b, c, d] * aux[a, b, c, d, e, f]

        for a, b, c, d in edges_source_target:
            source = sources[a - 1]
            target = targets[abs(c) - 1]
            L_out, U_out = eM.estima_M_complete(target, source)
            for e, f in vertices_target:
                model.addConstr(paux[a, b, c, d, e, f] >= L_out * aux[a, b, c, d, e, f])
                model.addConstr(paux[a, b, c, d, e, f] >= dist[a, b, c, d] - U_out * (1 - aux[a, b, c, d, e, f]))

        for a, b, c, d in edges_total:
            model.addConstr(p[a, b, c, d] == gp.quicksum(abs(e)*paux[a, b, c, d, e, f] for e, f in vertices_target))
        
    else:
        for a, b, c, d in edges_source:
            source = sources[a - 1]
            punto = neigh.Punto(barriers[c-1000][d])
            L_out, U_out = eM.estima_M_complete(source, punto)
            for g, h in vertices_target:
                model.addConstr(paux[a, b, c, d, a, 0, g, h] >= L_out * flow[a, b, c, d, a, 0, g, h])
                model.addConstr(paux[a, b, c, d, a, 0, g, h] >= dist[a, b, c, d] - U_out * (1 - flow[a, b, c, d, a, 0, g, h]))


        for a, b, c, d in edges_target:
            target = targets[abs(c) - 1]
            punto = neigh.Punto(barriers[a-1000][b])
            L_out, U_out = eM.estima_M_complete(target, punto)
            for e, f in vertices_source:
                model.addConstr(paux[a, b, c, d, e, f, c, 0] >= L_out * flow[a, b, c, d, e, f, c, 0])
                model.addConstr(paux[a, b, c, d, e, f, c, 0] >= dist[a, b, c, d] - U_out* (1 -flow[a, b, c, d, e, f, c, 0]))

        for a, b, c, d in edges_barrier:
            for e, f in vertices_source:
                for g, h in vertices_target:
                    # L_out, U_out = eM.estima_M_complete(neigh.Punto(np.array(barriers[a - 1000][b])), neigh.Punto(np.array(barriers[c - 1000][d])))
                    paux[a, b, c, d, e, f, g, h] = dist[a, b, c, d] * flow[a, b, c, d, e, f, g, h]
                    # model.addConstr(paux[a, b, c, d, e, f, g, h] >= L_out * flow[a, b, c, d, e, f, g, h])
                    # model.addConstr(paux[a, b, c, d, e, f, g, h] >= dist[a, b, c, d] - U_out * (1 - flow[a, b, c, d, e, f, g, h]))

        for a, b, c, d in edges_source_target:
            source = sources[a - 1]
            target = targets[abs(c) - 1]
            L_out, U_out = eM.estima_M_complete(target, source)
            model.addConstr(paux[a, b, c, d, a, 0, c, 0] >= L_out * flow[a, b, c, d, a, 0, c, 0])
            model.addConstr(paux[a, b, c, d, a, 0, c, 0] >= dist[a, b, c, d] - U_out * (1 - flow[a, b, c, d, a, 0, c, 0]))
        
        suma = 0

        for e, f in vertices_source:
            for g, h in vertices_target:
                model.addConstr(p[e, f, g, h] == paux.sum('*', '*', '*', '*', e, f, g, h))
    
    print('Spent time: ' + str(time.time() - start_product))
    
    print('\nSetting k-median constraints: ')
    start_median = time.time()

    # We only open k facilities
    model.addConstr(y.sum('*', '*') == k)

    # All the targets must be associated to one facility
    model.addConstrs(x.sum('*', '*', c, d) == 1 for c, d in vertices_target)

    # If a facility is not open in S, we can not launch capacity from this source.
    model.addConstrs(x[a, b, c, d] <= y[a, b] for a, b, c, d in x_index)

    print('Spent time: ' + str(time.time() - start_median))

    print('\nBinarying variables: ')
    start_binary = time.time()
    # Binarying an integer variable
    if single:
        for a, b, c, d in edges_total:
            model.addConstr(flow[a, b, c, d] == gp.quicksum((abs(e))*aux[a, b, c, d, e, f] for e, f in vertices_target))
    
    print('Spent time: ' + str(time.time() - start_binary))

    print('\nSetting flow constraints: ')
    start_flow = time.time()
    if single:
        for a, b in vertices_total:
            if (a, b) in vertices_source:
                model.addConstr(flow.sum(a, b, '*', '*') == x.sum(a, b, '*', '*'))
            elif (a, b) in vertices_barrier:
                model.addConstr(flow.sum(a, b, '*', '*') - flow.sum('*', '*', a, b) == 0)
            else:
                model.addConstr(flow.sum('*', '*', a, b) == 1)
    else:
        for e, f in vertices_source:
            for g, h in vertices_target:
                for a, b in vertices_total:
                    if (a, b) in vertices_source:
                        model.addConstr(flow.sum(e, f, '*', '*', e, f, g, h) == x[e, f, g, h])
                    elif (a, b) in vertices_barrier:
                        model.addConstr(flow.sum(a, b, '*', '*', e, f, g, h) - flow.sum('*', '*', a, b, e, f, g, h) == 0)
                    else:
                        model.addConstr(flow.sum('*', '*', g, h, e, f, g, h) == x[e, f, g, h])

    print('Spent time: ' + str(time.time() - start_flow))
# (1081, 0, 1083, 0, 3, 0, -6, 0),
    model.addConstr(flow[1081, 0, 1083, 1, 3, 0, -6, 0] + flow[1084, 1, 1083, 0, 3, 0, -6, 0] + flow[1084, 1, 1082, 1, 3, 0, -6, 0] + flow[1081, 0, 1083, 0, 3, 0, -6, 0] + flow[1081, 0, 1082, 1, 3, 0, -6, 0]<= 0)
    
    model.update()

    # Objective function
    objective = gp.quicksum(p[index] for index in p.keys()) + gp.quicksum(0.5*wL*flow[index] for index in flow.keys())
    model.setObjective(objective, GRB.MINIMIZE)

    # Model Setup Time
    second_time = time.time()
    time_elapsed = second_time - first_time

    model.update()

    # Model Params
    model.Params.Threads = 6
    model.Params.timeLimit = time_limit
    model.Params.NumericFocus = 1
    
    if lazy and single:
        model._flow = flow
        model._point = point
        model._alpha = alpha
        model._beta = beta
        model._gamma = gamma
        model.Params.LazyConstraints = 1

    if not (A4):
        model.Params.NonConvex = 2

    if lazy:
        model.optimize(elimcuts)
    elif init:
        model.optimize(callback)
    else:
        model.optimize()

    second_time = time.time()

    time_elapsed = second_time - first_time

    results = [len(sources), len(barriers), k, wL, lazy, A4, init, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]


    if init:
        try:
            results[-4] = time_h
            results[-2] = objval_h
            results[-3] = model._starttime
            results[-1] = model._startobjval
        except:
            print('No solution obtained by the heuristic')

    if model.Status == 3:
        model.computeIIS()
        model.write('infeasible_constraints.ilp')
        return results

    if model.SolCount == 0:
        return results

    # model.write('solution.sol')

    results[7] = model.getAttr('MIPGap')
    results[8] = model.getAttr('Runtime')
    results[9] = model.getAttr('NodeCount')
    results[10] = model.ObjVal

    if single:
        x_indices = [(index, x[index].X) for index in x.keys() if x[index].X > 0.5]
        # dist_indices = [(index, dist[index[0:4]].X) for index in flow.keys() if (flow[index].X > 0.5)]
        # p_indices = [(index, p[index].X) for index in p.keys() if flow[index].X > 0.5]
    else:

        x_indices = [(index, x[index].X) for index in x.keys() if x[index].X > 0.5]
        # dist_indices = [(index, dist[index[0:4]].X) for index in flow.keys() if flow[index].X > 0.5]
        # p_indices = [(index, p[index].X) for index in p.keys() if x[index].X > 0.5]

    x_indices = [index for index in x.keys() if x[index].X > 0.5]

    y_indices = [index for index in y.keys() if y[index].X > 0.5]

    flow_indices = [index for index in flow.keys() if flow[index].X > 0.5]

    if log:
        print(x_indices)
    
    # print(y_indices)
    
    print(flow_indices)

    if picture:
        # fig, ax = plt.subplots()

        # for b in barriers:
        #     ax.plot([b[0][0], b[1][0]], [b[0][1], b[1][1]], c='red')

        # for n in sources + targets:
        #     ax.add_artist(n.artist)

        points_vals = model.getAttr('x', point)
        print(points_vals)

        # for a, b in y_indices:
        #     ax.scatter(point[a, b, 0].X, point[a, b, 1].X, s = 30, c = 'black')

        segments = []

        if single:
            for a, b, c, d in flow_indices:
                if (a, b, c, d) in edges_source:
                    segments.append([point[a, b, 0].X, barriers[c-1000][d][0], point[a, b, 1].X, barriers[c-1000][d][1]])

                if (a, b, c, d) in edges_barrier:
                    segments.append([barriers[a-1000][b][0], barriers[c-1000][d][0], barriers[a-1000][b][1], barriers[c-1000][d][1]])

                if (a, b, c, d) in edges_target:
                    segments.append([barriers[a-1000][b][0], point[c, d, 0].X, barriers[a-1000][b][1], point[c, d, 1].X])

                if (a, b, c, d) in edges_source_target:
                    segments.append([point[a, b, 0].X, point[c, d, 0].X, point[a, b, 1].X, point[c, d, 1].X])
        else:
            for a, b, c, d, e, f, g, h in flow_indices:
                if (a, b, c, d) in edges_source:
                    segments.append([point[a, b, 0].X, barriers[c-1000][d][0], point[a, b, 1].X, barriers[c-1000][d][1]])

                if (a, b, c, d) in edges_barrier:
                    segments.append([barriers[a-1000][b][0], barriers[c-1000][d][0], barriers[a-1000][b][1], barriers[c-1000][d][1]])

                if (a, b, c, d) in edges_target:
                    segments.append([barriers[a-1000][b][0], point[c, d, 0].X, barriers[a-1000][b][1], point[c, d, 1].X])

                if (a, b, c, d) in edges_source_target:
                    segments.append([point[a, b, 0].X, point[c, d, 0].X, point[a, b, 1].X, point[c, d, 1].X])

        
        print(segments)
        # for segment in segments:
        #     ax.arrow(segment[0], segment[2], segment[1] - segment[0], segment[3] - segment[2], width=0.1,
        #              head_width=1, length_includes_head=True, color='black')


        # plt.axis([bounds_min, bounds_max, bounds_min, bounds_max])

        # ax.set_aspect('equal')
        # plt.show()

    return results
