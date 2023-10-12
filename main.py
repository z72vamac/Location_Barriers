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

from sflpn_b import sflpn_b
from h_kmedian_n import h_kmedian_n

# from HTSPS_without_prepro import HTSPS_without_prepro
# 50-3

# segments = np.genfromtxt('./instancias/segmentos50-3.csv', delimiter = ',')

# barriers = []
# for lista in segments:
#     barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

# bolas = np.genfromtxt('./instancias/bolas50-3.csv', delimiter = ',')
# sources = [neigh.Circle(center = [centro1, centro2], radii = radio, col = 'blue') for centro1, centro2, radio in bolas]

# targets = [neigh.Circle(center = [centro1, centro2], radii = radio, col = 'red') for centro1, centro2, radio in bolas]
#
# k = 3
#
# resultados = h_kmedian_n(barriers, sources, targets, k, time_limit=3600, prepro=True, log=False, picture=True)

# from h_kmedian_n_moresubindex import h_kmedian_n
# from h_kmedian_n import h_kmedian_n
#
## EXAMPLE NON-VISIBLE

barrier1 = [[0, 90], [30, 60]]
barrier2 = [[10, 50], [40, 50]]
barrier3 = [[0, 30], [10, 40]]
barrier4 = [[10, 30], [30, 5]]
barrier5 = [[40, 10], [70, 40]]
barrier6 = [[60, 20], [100, 10]]
barrier7 = [[30, 70], [70, 95]]
barrier8 = [[70, 90], [60, 50]]
barrier9 = [[70, 80], [90, 60]]
barrier10 = [[74, 33], [98, 60]]

barriers = [barrier1, barrier2, barrier3, barrier4, barrier5, barrier6, barrier7, barrier8, barrier9, barrier10]
# barriers = [barrier3, barrier5, barrier6, barrier7, barrier8]

N1s = neigh.Circle(center=[70, 55], radii=4, col = 'green')
N2s = neigh.Circle(center=[50, 70], radii=8, col = 'green')
N3s = neigh.Circle(center=[30, 35], radii=10, col = 'green')

sources = [N1s, N2s, N3s]

N1t = neigh.Circle(center=[10, 15], radii=6, col = 'blue')
N2t = neigh.Circle(center=[65, 10], radii=7, col = 'blue')
N3t = neigh.Circle(center=[10, 65], radii=5, col = 'blue')
N4t = neigh.Circle(center=[90, 35], radii=6, col = 'blue')
N5t = neigh.Circle(center=[90, 85], radii=6, col = 'blue')
N6t = neigh.Circle(center=[30, 90], radii=10, col = 'blue')

targets = [N1t, N2t, N3t, N4t, N5t, N6t]

k = 2

endurance = 1000

wE = 1

wL = 0

# 279.88
resultados = h_kmedian_n(barriers, sources=sources, targets=targets, k=k, wL=wL, single=False, lazy=False, A4=False, time_limit=300, picture=True, init = False)

# print(resultados)

# np.random.seed(5)

# instances = range(5)
# n_Ns = [10, 20, 30, 50, 80]
# # n_Ns = [25, 30]

# ks = []

# for nn in n_Ns:
#     ks.append(k)

# time_limit = 3600

# wE = 1

# # wLs = [50, 100]
# wLs = [50]

# lazy = [False]
# # A4 = [True, False]
# A4 = [True, False]

# inits = [False, True]
# # inits = [True]

# singles = [False, True]
# time_limit = 3600

# instance = 0

# nn = 10

# bolas = np.genfromtxt('./instancias/bolas' + str(nn) + '-' + str(instance) + '.csv', delimiter=',')

# N = [neigh.Circle(center=[centro1, centro2], radii=radio, col = 'blue') for centro1, centro2, radio in bolas]

# segments = np.genfromtxt('./instancias/segmentos' + str(nn) + '-' + str(instance) + '.csv', delimiter=',')

# barriers = []

# for lista in segments:
#     barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

# nB = len(barriers)

# sublist = np.random.choice(nB, int(np.floor(0.5*nB)))
# barriers1 = [barriers[b] for b in sublist]

# resultados = h_kmedian_n(barriers1, sources=N, targets=N, k=k, single=True, wL=50, lazy=False, A4=False, init=False, picture = False,
#                             time_limit=3600) 