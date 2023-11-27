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
from h_kmedian_new_casestudy import h_kmedian_n


blocks = [9]
barriers = []
sources = []
targets = []

for i in blocks:

    segments = np.genfromtxt('./instances_casestudy/barriers/barriers{0}.csv'.format(i), delimiter = ',')

    for lista in segments:
        barriers.append([[lista[0], lista[1]], [lista[2], lista[3]]])

    neighbourhoods = np.genfromtxt('./instances_casestudy/circles/sources/circles{0}.csv'.format(i), delimiter = ',')

    for x_coord, y_coord, radii in neighbourhoods:
        sources.append(neigh.Circle(center = [x_coord, y_coord], radii = radii, col='blue'))

    neighbourhoods = np.genfromtxt('./instances_casestudy/circles/targets/circles{0}.csv'.format(i), delimiter = ',')

    for x_coord, y_coord, radii in neighbourhoods:
        targets.append(neigh.Circle(center = [x_coord, y_coord], radii = radii, col="green"))

# bolas = [[59.5, 57.5, 7.5], [(57+65.5)/2, 200-(92+100.5)/2, 4.25], [(77+81)/2, 200-82, 2], [(98+103)/2, 200-(85+90)/2, 2.5], [(117+124)/2, 200-(100+107)/2, 3.5]]

# N = [neigh.Circle(center = [centro1, centro2], radii = radio) for centro1, centro2, radio in bolas]

for k in range(1, 5):
    resultados = h_kmedian_n(barriers, sources=sources, targets=targets, k=3, wL=0, single=False, lazy=False, A4=False, time_limit=3600, picture=True, init = False)
