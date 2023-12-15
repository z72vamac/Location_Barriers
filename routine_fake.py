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
from h_kmedian_new import h_kmedian_n
from neighborhood import Circle
import sys

# Parameters of the model
# instances = range(5)


instance = sys.argv[1]

print("Solving problems of instance:" + str(instance))

n_Ns = [10, 20, 30, 50, 70]
# n_Ns = [25, 30]

ks = []
percs_k = []

for nn in n_Ns:
    k = [1, math.floor(0.1*nn), math.floor(0.25*nn)]
    ks.append(k)
    
    perc_k = [1, '10%', '25%']
    percs_k.append(perc_k)

time_limit = 3600

wE = 1

wLs = [0, 50]

lazy = [False]
percs = [0.1, 0.2, 0.5, 1]

inits = [False, True]

singles = [True, False]

time_limit = 3600

start = False


if start:
    dataframe = pd.read_csv('./resultados/parameters_' + str(instance) + '.csv').iloc[:, 1:]
    num_rows = dataframe.shape[0]
else:
    num_rows = 0
    dataframe = pd.DataFrame(columns=['instance', 'n_N', 'perc_B', 'k', 'single', 'wL', 'lazy', 'A4', 'init'])

counter = 1


for nn, perc_k in zip(n_Ns, percs_k):
    for single in singles:
            for k in perc_k:
                for wL in wLs:
                    for perc in percs:
                        for l in lazy:
                            for init in inits:
                                circles = np.genfromtxt('./instances_random/circles/circles' + str(nn) + '-' + str(instance) + '.csv', delimiter=',')
                                neighbourhoods = [Circle(center=[centro1, centro2], radii=radio, col='blue') for centro1, centro2, radio in circles]

                                segments = np.genfromtxt('./instances_random/barriers/barriers' + str(nn) + '-' + str(instance) + '.csv', delimiter=',')
                                barriers = []

                                for lista in segments:
                                    barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

                                nB = len(barriers)
                                
                                np.random.seed(5)
                                # sublist = np.random.choice(nB, int(np.floor(perc * nB)))
                                # barriers1 = [barriers[b] for b in sublist]

                                if counter > num_rows:

                                    print('\n\nSolving hampered k-median')
                                    print('Instance: ' + str(instance))
                                    print('Number of neighbourhoods: ' + str(nn))
                                    print('k: ' + str(k))
                                    print('Lazy mode: ' + str(l))
                                    print('Init: ' + str(init))
                                    print('Percentage of barriers: ' + str(perc) + '%\n\n')


                                    if perc < 1:
                                        A4 = False
                                    else:
                                        A4 = True

                                    # resultados = h_kmedian_n(barriers1, sources=neighbourhoods, targets=neighbourhoods, k=k, single=single, wL=wL, lazy=l, A4=A4, init=init, time_limit=time_limit)
                                    perc2 = str(perc*100) + '%' 
                                    serie = pd.Series([instance, nn, perc2, k, single, wL, l, A4, init], index=dataframe.columns)

                                    dataframe = dataframe._append(serie, ignore_index=True)
                                    dataframe.to_csv('./resultados/parameters_' + str(instance) + '.csv')

                                counter += 1