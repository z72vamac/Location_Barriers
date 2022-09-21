# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 00:16:23 2019

@author: carlo
"""
from gurobipy import *
import numpy as np
from neighborhood import *
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import auxiliar_functions as af


# def estima_BigM(data):
#
#     m = len(data)
#     BigM = 0
#
#     for i in range(m):
#         for j in range(m):
#             if i != j:
#                 comp1 = data[i]
#                 comp2 = data[j]
#
#                 if type(comp1) is Poligono or Poligonal:
#                     if type(comp2) is Poligono:
#                         maximo = max([np.linalg.norm(v-w) for v in comp1.V
#                                                           for w in comp2.V])
#
#                     if type(comp2) is Elipse:
#                         maximo = comp2.radio + max([np.linalg.norm(v-comp2.centro) for v in comp1.V])
#
#                 if type(comp1) is Poligonal:
#                     if type(comp2) is Poligono or Poligonal:
#                         maximo = max([np.linalg.norm(w-v) for w in comp2.V
#                                                           for v in comp1.V])
#
#                     if type(comp2) is Elipse:
#                         caso1 = comp2.radio + np.linalg.norm(comp1.P-comp2.centro)
#                         caso2 = comp2.radio + np.linalg.norm(comp1.Q-comp2.centro)
#                         maximo = max(caso1, caso2)
#
#                 if type(comp1) is Elipse:
#                     if type(comp2) is Poligono or Poligonal:
#                         maximo = comp1.radio + max([np.linalg.norm(comp1.centro-w) for w in comp2.V])
#
#                     if type(comp2) is Elipse:
#                         maximo = comp1.radio + np.linalg.norm(comp1.centro-comp2.centro) + comp2.radio
#
#                 if maximo >= BigM:
#                     BigM = maximo
#
#     return BigM


def estima_BigM_local(comp1, comp2):
        maximo = 0
        if type(comp1) is Poligono or type(comp1) is Poligonal:
            if type(comp2) is Poligono or type(comp2) is Poligonal:
                maximo = max([np.linalg.norm(v-w) for v in comp1.V
                                                  for w in comp2.V])

            if type(comp2) is Elipse:
                maximo = comp2.radio + max([np.linalg.norm(v-comp2.centro) for v in comp1.V])

        if type(comp1) is Elipse:
            if type(comp2) is Poligono or type(comp2) is Poligonal:
                maximo = comp1.radio + max([np.linalg.norm(comp1.centro-w) for w in comp2.V])

            if type(comp2) is Elipse:
                maximo = comp1.radio + np.linalg.norm(comp1.centro-comp2.centro) + comp2.radio

        return maximo

def estima_SmallM_local(comp1, comp2):
        if type(comp1) is Poligono or type(comp1) is Poligonal:
            if type(comp2) is Poligono or type(comp2) is Poligonal:
                minimo = min([np.linalg.norm(v-w) for v in comp1.V
                                                  for w in comp2.V])

            if type(comp2) is Elipse:
                minimo = - comp2.radio + min([np.linalg.norm(v-comp2.centro) for v in comp1.V])

        if type(comp1) is Elipse:
            if type(comp2) is Poligono or type(comp2) is Poligonal:
                minimo = -comp1.radio + min([np.linalg.norm(comp1.centro-w) for w in comp2.V])

            if type(comp2) is Elipse:
                minimo = -comp1.radio + np.linalg.norm(comp1.centro-comp2.centro) - comp2.radio

        return minimo

def estima_max_inside(comp):
        maximo = 0
        if type(comp) is Poligono:
            maximo = max([np.linalg.norm(v-w) for v in comp.V for w in comp.V])

        if type(comp) is Elipse:
            maximo = 2*comp.radio

        if type(comp) is Poligonal:
            maximo = comp.alpha * comp.longitud

        return maximo

def estima_M_alpha1(entorno, punto1, punto2):
    if type(entorno) is Circle:

        theta = np.linspace(0, 2*np.pi, 100)

        centro = entorno.center
        radio = entorno.radii

        x = centro[0] + radio*np.cos(theta)
        y = centro[1] + radio*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant([x[i], y[i]], punto1, punto2) for i in range(100)]

        m = -100*abs(min(determinantes))
        M = 100*abs(max(determinantes))

        return m, M

def estima_M_alpha2(punto1, entorno, punto2):
    if type(entorno) is Circle:

        theta = np.linspace(0, 2*np.pi, 100)

        centro = entorno.center
        radio = entorno.radii

        x = centro[0] + radio*np.cos(theta)
        y = centro[1] + radio*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant(punto1, [x[i], y[i]], punto2) for i in range(100)]

        m = -100*abs(min(determinantes))
        M = 100*abs(max(determinantes))

        return m, M

def estima_M_alpha3(punto1, punto2, entorno):
    if type(entorno) is Circle:

        theta = np.linspace(0, 2*np.pi, 100)

        centro = entorno.center
        radio = entorno.radii

        x = centro[0] + radio*np.cos(theta)
        y = centro[1] + radio*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant(punto1, punto2, [x[i], y[i]]) for i in range(100)]

        m = -100*abs(min(determinantes))
        M = 100*abs(max(determinantes))

        return m, M

def estima_M_alpha4(punto1, entorno1, entorno2):
    if type(entorno1) is Circle and type(entorno2) is Circle:

        theta = np.linspace(0, 2*np.pi, 100)

        centro1 = entorno1.center
        radio1 = entorno1.radii

        x1 = centro1[0] + radio1*np.cos(theta)
        y1 = centro1[1] + radio1*np.sin(theta)

        centro2 = entorno2.center
        radio2 = entorno2.radii

        x2 = centro2[0] + radio2*np.cos(theta)
        y2 = centro2[1] + radio2*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant(punto1, [x1[i], y1[i]], [x2[j], y2[j]]) for i in range(100) for j in range(100)]

        m = -100*abs(min(determinantes))
        M = 100*abs(max(determinantes))

        return m, M
