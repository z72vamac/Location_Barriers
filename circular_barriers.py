import time
import itertools

import gurobipy as gp
from gurobipy import GRB
from matplotlib.patches import Circle

import auxiliar_functions as af
from data import *
import neighborhood as neigh


neighborhood1 = neigh.Circle(center=[4, 4], radii=1, col='green')
neighborhood2 = neigh.Circle(center=[15, 11], radii=1, col='green')

barrier = neigh.Circle(center = [9, 8], radii = 4, col = 'red')

model = gp.Model('Model: C-KMedian-N')

mid1 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'mid1')
mid2 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'mid2')

pn1 = model.addVars(2, vtype = GRB.CONTINUOUS, name = "pn1")
pn2 = model.addVars(2, vtype = GRB.CONTINUOUS, name = "pn2")

pb1 = model.addVars(2, vtype = GRB.CONTINUOUS, name = "pb1")
pb2 = model.addVars(2, vtype = GRB.CONTINUOUS, name = "pb2")

dn1 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dn1')
difn1 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'difn1')
dn2 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dn2')
difn2 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'difn2')
db12 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'db12')
difb12 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'difb12')

dfrontn1 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dfrontn1')
diffrontn1 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'diffrontn1')

dfrontb1 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dfrontb1')
diffrontb1 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'diffrontb1')
dfrontmidb1 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dfrontmidb1')
diffrontmidb1 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'diffrontmidb1')


dfrontb2 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dfrontb2')
diffrontb2 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'diffrontb2')
dfrontmidb2 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dfrontmidb2')
diffrontmidb2 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'diffrontmidb2')

dfrontn2 = model.addVar(vtype = GRB.CONTINUOUS, lb = 0.0, name = 'dfrontn2')
diffrontn2 = model.addVars(2, vtype = GRB.CONTINUOUS, name = 'diffrontn2')

model.update()

# Creacion de los puntos medios
model.addConstr(2*mid1[0] == pn1[0] + barrier.center[0])
model.addConstr(2*mid1[1] == pn1[1] + barrier.center[1])

model.addConstr(2*mid2[0] == pn2[0] + barrier.center[0])
model.addConstr(2*mid2[1] == pn2[1] + barrier.center[1])

# Pn1 esta en la bola 1
model.addConstr(diffrontn1[0] >=   pn1[0] - neighborhood1.center[0])
model.addConstr(diffrontn1[0] >= - pn1[0] + neighborhood1.center[0])

model.addConstr(diffrontn1[1] >=   pn1[1] - neighborhood1.center[1])
model.addConstr(diffrontn1[1] >= - pn1[1] + neighborhood1.center[1])

model.addConstr(diffrontn1[0]*diffrontn1[0] + diffrontn1[1]*diffrontn1[1] <= dfrontn1*dfrontn1)

# Pn2 esta en la bola2
model.addConstr(diffrontn2[0] >=   pn2[0] - neighborhood2.center[0])
model.addConstr(diffrontn2[0] >= - pn2[0] + neighborhood2.center[0])

model.addConstr(diffrontn2[1] >=   pn2[1] - neighborhood2.center[1])
model.addConstr(diffrontn2[1] >= - pn2[1] + neighborhood2.center[1])

model.addConstr(diffrontn2[0]*diffrontn2[0] + diffrontn2[1]*diffrontn2[1] <= dfrontn2*dfrontn2)

# Pb1 esta en la bola barrera
model.addConstr(diffrontb1[0] >=   pb1[0] - barrier.center[0])
model.addConstr(diffrontb1[0] >= - pb1[0] + barrier.center[0])

model.addConstr(diffrontb1[1] >=   pb1[1] - barrier.center[1])
model.addConstr(diffrontb1[1] >= - pb1[1] + barrier.center[1])

model.addConstr(diffrontb1[0]*diffrontb1[0] + diffrontb1[1]*diffrontb1[1] <= dfrontb1*dfrontb1)

model.addConstr(dfrontb1 <= barrier.radii)

# Pb1 esta en la bola centrada en el punto medio
model.addConstr(diffrontmidb1[0] >=   pb1[0] - mid1[0])
model.addConstr(diffrontmidb1[0] >= - pb1[0] + mid1[0])

model.addConstr(diffrontmidb1[1] >=   pb1[1] - mid1[1])
model.addConstr(diffrontmidb1[1] >= - pb1[1] + mid1[1])

model.addConstr(diffrontmidb1[0]*diffrontmidb1[0] + diffrontmidb1[1]*diffrontmidb1[1] <= dfrontmidb1*dfrontmidb1)

# Pb2 esta en la bola barrera
model.addConstr(diffrontb2[0] >=   pb2[0] - barrier.center[0])
model.addConstr(diffrontb2[0] >= - pb2[0] + barrier.center[0])

model.addConstr(diffrontb2[1] >=   pb2[1] - barrier.center[1])
model.addConstr(diffrontb2[1] >= - pb2[1] + barrier.center[1])

model.addConstr(diffrontb2[0]*diffrontb2[0] + diffrontb2[1]*diffrontb2[1] <= dfrontb2*dfrontb2)

model.addConstr(dfrontb1 <= barrier.radii)

# Pb2 esta en la bola centrada en el punto medio
model.addConstr(diffrontmidb2[0] >=   pb2[0] - mid2[0])
model.addConstr(diffrontmidb2[0] >= - pb2[0] + mid2[0])

model.addConstr(diffrontmidb2[1] >=   pb2[1] - mid2[1])
model.addConstr(diffrontmidb2[1] >= - pb2[1] + mid2[1])

model.addConstr(diffrontmidb2[0]*diffrontmidb2[0] + diffrontmidb2[1]*diffrontmidb2[1] <= dfrontmidb2*dfrontmidb2)

# Distancia entre pn1 y pb1
model.addConstr(difn1[0] >=   pn1[0] - pb1[0])
model.addConstr(difn1[0] >= - pn1[0] + pb1[0])

model.addConstr(difn1[1] >=   pn1[1] - pb1[1])
model.addConstr(difn1[1] >= - pn1[1] + pb1[1])

model.addConstr(difn1[0]*difn1[0] + difn1[1]*difn1[1] <= dn1*dn1)

# Distancia entre pb1 y pb2
model.addConstr(difb12[0] >=   pb1[0] - pb2[0])
model.addConstr(difb12[0] >= - pb1[0] + pb2[0])

model.addConstr(difb12[1] >=   pb1[1] - pb2[1])
model.addConstr(difb12[1] >= - pb1[1] + pb2[1])

model.addConstr(difb12[0]*difb12[0] + difb12[1]*difb12[1] <= db12*db12)

# Distancia entre pn2 y pb2
model.addConstr(difn2[0] >=   pn2[0] - pb2[0])
model.addConstr(difn2[0] >= - pn2[0] + pb2[0])

model.addConstr(difn2[1] >=   pn2[1] - pb2[1])
model.addConstr(difn2[1] >= - pn2[1] + pb2[1])

model.addConstr(difn2[0]*difn2[0] + difn2[1]*difn2[1] <= dn2*dn2)













