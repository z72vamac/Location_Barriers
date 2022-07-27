import gurobipy as gp
from gurobipy import GRB
import numpy as np
from itertools import product, permutations, chain
import random
import matplotlib.pyplot as plt
# from matplotlib.patches import Circle, Polygon
from matplotlib.collections import PatchCollection
from data import *
import neighborhood as neigh
import copy
import estimacion_M as eM
import auxiliar_functions as af
import networkx as nx
import math
import numpy as np
import pandas as pd
from h_kmedian_n import h_kmedian_n
from neighborhood import Circle

# Parameters of the model
instances = range(5)
n_Ns = [10, 20, 30, 50, 80, 100]
ks = []

for nn in n_Ns:
    k = [1, math.floor(0.1*nn), math.floor(0.25*nn)]
    ks.append(k)

time_limit = 3600

wE = 1

wL = 50

lazy = [False, True]
A4 = [True, False]

init = False

dataframe = pd.DataFrame(columns=['Instance', 'n_N', 'n_B', 'k', 'wL', 'Lazy', 'A4', 'Gap', 'Runtime', 'NodeCount', 'ObjVal', 'Runtime_h', 'ObjVal_h'])

for instance in instances:
    for nn, k_list in zip(n_Ns, ks):
        for k in k_list:
            for a in A4:
                for l in lazy:


                    bolas = np.genfromtxt('./instancias/bolas' + str(nn) + '-' + str(instance) + '.csv', delimiter=',')

                    N = [Circle(center=[centro1, centro2], radii=radio, col = 'blue') for centro1, centro2, radio in bolas]

                    segments = np.genfromtxt('./instancias/segmentos' + str(nn) + '-' + str(instance) + '.csv', delimiter=',')

                    barriers = []

                    nB = len(barriers)

                    for lista in segments:
                        barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

                    if not(a):
                        print('\n\n Solving hampered k-median')
                        print('Instance: ' + str(instance))
                        print('Number of neighbourhoods: ' + str(nn))
                        print('Instances verifying A4: ' + str(a))
                        print('Lazy mode: ' + str(l))
                        print('Percentage of barriers: 10 %\n\n')

                        barriers = np.random.choice(barriers, size = np.floor(0.1*nB), replace=False)

                        resultados = h_kmedian_n(barriers, sources=N, targets=N, k=k, wL=wL, lazy=l, A4=a,
                                                 time_limit=time_limit)

                        serie = pd.Series([instance] + resultados, index=dataframe.columns)

                        dataframe = dataframe.append(serie, ignore_index=True)
                        dataframe.to_csv('./resultados/results.csv')

                        print('\n\n Solving hampered k-median')
                        print('Instance: ' + str(instance))
                        print('Number of neighbourhoods: ' + str(nn))
                        print('Instances verifying A4: ' + str(a))
                        print('Lazy mode: ' + str(l))
                        print('Percentage of barriers: 20 %\n\n')

                        barriers = np.random.choice(barriers, size = np.floor(0.2*nB), replace=False)

                        resultados = h_kmedian_n(barriers, sources=N, targets=N, k=k, wL=wL, lazy=l, A4=a,
                                                 time_limit=time_limit)

                        serie = pd.Series([instance] + resultados, index=dataframe.columns)

                        dataframe = dataframe.append(serie, ignore_index=True)
                        dataframe.to_csv('./resultados/results.csv')

                        print('\n\n Solving hampered k-median')
                        print('Instance: ' + str(instance))
                        print('Number of neighbourhoods: ' + str(nn))
                        print('Instances verifying A4: ' + str(a))
                        print('Lazy mode: ' + str(l))
                        print('Percentage of barriers: 50 %\n\n')

                        barriers = np.random.choice(barriers, size = np.floor(0.5*nB), replace=False)

                        resultados = h_kmedian_n(barriers, sources=N, targets=N, k=k, wL=wL, lazy=l, A4=a,
                                                 time_limit=time_limit)

                        serie = pd.Series([instance] + resultados, index=dataframe.columns)

                        dataframe = dataframe.append(serie, ignore_index=True)
                        dataframe.to_csv('./resultados/results.csv')

                    else:
                        print('Percentage of barriers: 100 %\n\n')

                        resultados = h_kmedian_n(barriers, sources=N, targets=N, k=k, wL=wL, lazy=l, A4=a,
                                                 time_limit=time_limit)

                        serie = pd.Series([instance] + resultados, index=dataframe.columns)

                        dataframe = dataframe.append(serie, ignore_index=True)
                        dataframe.to_csv('./resultados/results.csv')