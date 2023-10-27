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

def preproM(m, M, factor=2):
    if m < 0 and M <= 0:
        m *= factor
        M /= factor
    elif m < 0 and M > 0:
        m *= factor
        M *= factor
    else:
        m /= factor
        M *= factor

    return m, M

n_iter = 100

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

def estima_M_alpha1(entorno, punto1, punto2, n_iter = 20):
    if type(entorno) is Circle:

        theta = np.linspace(0, 2*np.pi, n_iter)

        centro = entorno.center
        radio = entorno.radii

        x = centro[0] + radio*np.cos(theta)
        y = centro[1] + radio*np.sin(theta)

        determinantes = [af.determinant([x[i], y[i]], punto1, punto2) for i in range(n_iter)]

        m = min(determinantes)
        M = max(determinantes)

        m, M = preproM(m, M)

        return m, M

def estima_M_alpha2(punto1, entorno, punto2, n_iter = 20):
    if type(entorno) is Circle:

        theta = np.linspace(0, 2*np.pi, n_iter)

        centro = entorno.center
        radio = entorno.radii

        x = centro[0] + radio*np.cos(theta)
        y = centro[1] + radio*np.sin(theta)

        determinantes = [af.determinant(punto1, [x[i], y[i]], punto2) for i in range(n_iter)]

        m = min(determinantes)
        M = max(determinantes)

        m, M = preproM(m, M)

        return m, M

def estima_M_alpha3(punto1, punto2, entorno, n_iter = 20):
    if type(entorno) is Circle:

        theta = np.linspace(0, 2*np.pi, n_iter)

        centro = entorno.center
        radio = entorno.radii

        x = centro[0] + radio*np.cos(theta)
        y = centro[1] + radio*np.sin(theta)

        # print([x[0], y[0]])

        determinantes = [af.determinant(punto1, punto2, [x[i], y[i]]) for i in range(n_iter)]

        m = min(determinantes)
        M = max(determinantes)

        m, M = preproM(m, M)

        return m, M

def estima_M_alpha4(punto1, entorno1, entorno2, n_iter=20):
    if type(entorno1) is Circle and type(entorno2) is Circle:

        theta = np.linspace(0, 2*np.pi, n_iter)

        centro1 = entorno1.center
        radio1 = entorno1.radii

        x1 = centro1[0] + radio1*np.cos(theta)
        y1 = centro1[1] + radio1*np.sin(theta)

        centro2 = entorno2.center
        radio2 = entorno2.radii

        x2 = centro2[0] + radio2*np.cos(theta)
        y2 = centro2[1] + radio2*np.sin(theta)

        determinantes = [af.determinant(punto1, [x1[i], y1[i]], [x2[j], y2[j]]) for i in range(n_iter) for j in range(n_iter)]

        m = min(determinantes)
        M = max(determinantes)

        m, M = preproM(m, M)

        return m, M

def estima_M_complete(ent1, ent2, n_iter = 10):

    distancias = []

    if type(ent1) is Elipse and type(ent2) is Elipse:
        theta = np.linspace(0, 2*np.pi, n_iter)

        centro1 = ent1.center
        radii11 = ent1.width
        radii12 = ent1.height

        x1 = centro1[0] + radii11*np.cos(theta)
        y1 = centro1[1] + radii12*np.sin(theta)

        centro2 = ent2.center
        radii21 = ent2.width
        radii22 = ent2.height

        x2 = centro2[0] + radii21*np.cos(theta)
        y2 = centro2[1] + radii22*np.sin(theta)

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array([x2[j], y2[j]])) for i in range(n_iter) for j in range(n_iter)]

    elif type(ent1) is Elipse and type(ent2) is Circle:
        theta = np.linspace(0, 2*np.pi, n_iter)

        centro1 = ent1.center
        radii11 = ent1.width
        radii12 = ent1.height

        x1 = centro1[0] + radii11*np.cos(theta)
        y1 = centro1[1] + radii12*np.sin(theta)

        centro2 = ent2.center
        radii2 = ent2.radii

        x2 = centro2[0] + radii2*np.cos(theta)
        y2 = centro2[1] + radii2*np.sin(theta)

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array([x2[j], y2[j]])) for i in range(n_iter) for j in range(n_iter)]

    elif type(ent1) is Elipse and type(ent2) is Punto:
        theta = np.linspace(0, 2*np.pi, n_iter)

        centro = ent1.center
        radii11 = ent1.width
        radii12 = ent1.height

        x = centro[0] + radii11*np.cos(theta)
        y = centro[1] + radii12*np.sin(theta)

        distancias = [np.linalg.norm(np.array([x[i], y[i]]) - np.array(ent2.V)) for i in range(n_iter)]

    elif type(ent1) is Elipse and (type(ent2) is Poligono or type(ent2) is Poligonal):

        theta = np.linspace(0, 2*np.pi, n_iter)

        centro1 = ent1.center
        radii11 = ent1.width
        radii12 = ent1.height

        x1 = centro1[0] + radii11*np.cos(theta)
        y1 = centro1[1] + radii12*np.sin(theta)

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array(v)) for i in range(n_iter) for v in ent2.V]

    elif type(ent1) is Circle and type(ent2) is Elipse:
        return estima_M_complete(ent2, ent1)

    elif type(ent1) is Circle and type(ent2) is Circle:
        theta = np.linspace(0, 2*np.pi, n_iter)

        centro1 = ent1.center
        radii1 = ent1.radii

        x1 = centro1[0] + radii1*np.cos(theta)
        y1 = centro1[1] + radii1*np.sin(theta)

        centro2 = ent2.center
        radii2 = ent2.radii

        x2 = centro2[0] + radii2*np.cos(theta)
        y2 = centro2[1] + radii2*np.sin(theta)

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array([x2[j], y2[j]])) for i in range(n_iter) for j in range(n_iter)]

    elif type(ent1) is Circle and type(ent2) is Punto:
        theta = np.linspace(0, 2*np.pi, n_iter)

        centro = ent1.center
        radii1 = ent1.radii

        x = centro[0] + radii1*np.cos(theta)
        y = centro[1] + radii1*np.sin(theta)

        distancias = [np.linalg.norm(np.array(ent2.V) - np.array([x[i], y[i]])) for i in range(n_iter)]

    elif type(ent1) is Circle and (type(ent2) is Poligono or type(ent2) is Poligonal):
        theta = np.linspace(0, 2*np.pi, n_iter)

        centro1 = ent1.center
        radii1 = ent1.radii

        x1 = centro1[0] + radii1*np.cos(theta)
        y1 = centro1[1] + radii1*np.sin(theta)

        distancias = [np.linalg.norm(np.array([x1[i], y1[i]]) - np.array(v)) for i in range(n_iter) for v in ent2.V]

    elif type(ent1) is Punto and type(ent2) is Elipse:
        return estima_M_complete(ent2, ent1)

    elif type(ent1) is Punto and type(ent2) is Circle:
        return estima_M_complete(ent2, ent1)

    elif type(ent1) is Punto and type(ent2) is Punto:
        distancias = [np.linalg.norm(np.array(ent1.V) - np.array(ent2.V))]

        return min(distancias), max(distancias)

    elif type(ent1) is Punto and (type(ent2) is Poligono or type(ent2) is Poligonal):
        distancias = [np.linalg.norm(np.array(ent1.V) - np.array(v)) for v in ent2.V]
    
    elif (type(ent1) is Poligono or type(ent1) is Poligonal) and type(ent2) is Elipse:
        return estima_M_complete(ent2, ent1)

    elif (type(ent1) is Poligono or type(ent1) is Poligonal) and type(ent2) is Circle:
        return estima_M_complete(ent2, ent1)

    elif (type(ent1) is Poligono or type(ent1) is Poligonal) and type(ent2) is Punto:
        return estima_M_complete(ent2, ent1)

    elif (type(ent1) is Poligono or type(ent1) is Poligonal) and (type(ent2) is Poligono or type(ent2) is Poligonal):
        distancias = [np.linalg.norm(np.array(v) - np.array(w)) for v in ent1.V for w in ent2.V]

    else:
        "Llego aqui?"
    
    m = min(distancias)
    M = max(distancias)

    m, M = preproM(m, M)

    return m, M