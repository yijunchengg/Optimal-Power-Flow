#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:59:07 2022

@author: Yijun Cheng

"""

import numpy as np
import pandas as pd
import cmath
import time
from gurobipy import Model
from gurobipy import GRB
from gurobipy import quicksum


# data
branch_data = [
    [1,2,0.0922,0.047],
    [2,3,0.493,0.2511],
    [3,4,0.366,0.1864],
    ]
bus_data = [
    [0,0],
    [0.1,0.06],
    [0.09,0.04],
    [0.12,0.08],
    ]


branch = pd.DataFrame(columns=('bus1', 'bus2', 'resistance', 'reactance'),
                      data=branch_data)

load = pd.DataFrame(columns=('activepower', 'reactivepower'),
                    data=bus_data[:4]
                    )

n_b = len(branch)
n_l = len(load)

# Admittance matrix
Y = pd.DataFrame(data=np.zeros((n_l, n_l)))
for i in range(len(branch)):
    Y.iloc[branch.iloc[i, 0]-1, branch.iloc[i, 1]-1]=-1/complex(
        branch.iloc[i, 2], branch.iloc[i, 3])
    Y.iloc[branch.iloc[i, 1]-1, branch.iloc[i, 0]-1]=-1/complex(
        branch.iloc[i, 2], branch.iloc[i, 3])

for i in range(len(Y)):
    Y.iloc[i, i] = -sum(Y.loc[i])


start_time = time.time()
# load flow problem
# Create a new model
m = Model()

# Nonconvex problem
m.params.NonConvex = 2
# m.Params.TIME_LIMIT = 10

m.setParam('OutputFlag', 1)

# Create variables
# bus 1 is reference bus, slack bus
# voltage magnitude v and phase angle delta are knownas as 1 and 0
# active p and reactive q are unknown
s = m.addVars(2, lb=0, ub=4, name='substation')

# bus 3 is generator bus
# minimize active p, voltage magnitude should be known
# reactive q and phase angle delta are unknown
g_p = m.addVar(lb=0, ub=0.3, name='active power of generator')
# g_v = m.addVars(n_g, lb=0.95, ub=1.042, name='voltage magnitude of generator')
g_q = m.addVar(lb=0, ub=0.25, name='reactive power of generator')
# g_a = m.addVars(n_g, lb=-3.14, ub=3.14, name='phase angle of generator')


# load bus
# active p and reactive q are known
# voltage magnitude v and phase angle delta are unknown
l_v = m.addVars(n_l, lb=0.95, ub=1.042, name='voltage magnitude of load')
l_a = m.addVars(n_l, lb=-3.14, ub=3.14, name='phase angle of load')

m.addConstr(l_v[0]==1)
m.addConstr(l_a[0]==0)


# ************************************************
# Error
# addGenConstrCos() must be variables (Var), not linear expressions (LinExpr).
# introduce new variables t_cos
# Error
# Invalid argument to QuadExpr multiplication
# introduce new variable for solving by Gurobi
v_i_v_j = m.addVars(n_l, n_l, lb=0.95*0.95, ub=1.042*1.042,
                    name='v_i multiple v_j')

theta_ji= m.addVars(n_l, n_l, lb=-6.28, ub=6.28,
                      name='l_a[j]-l_a[i]')
theta_ij = m.addVars(n_l, n_l, lb=-6.28, ub=6.28,
                       name='l_a[i]-l_a[j]')

cos_theta_ji = m.addVars(n_l, n_l, lb=-6.28, ub=6.28,
                         name='cos(l_a[j]-l_a[i])')

theta = m.addVars(n_l, n_l, lb=-3.14, ub=3.14,
                  name='sum of theta')

cos_theta = m.addVars(n_l, n_l, lb=-1, ub=1,
                      name='cos of theta')
sin_theta = m.addVars(n_l, n_l, lb=-1, ub=1,
                      name='sin of theta')

for i in range(n_l):
    inject_active = 0
    inject_reactive = 0
    if i == 0:
        inject_active = s[0]
        inject_reactive = s[1]
    if i == 3:
        inject_active = g_p
        inject_reactive = g_q
    
    # sparse Y, if Y.iloc[i, j]!=0j
    for j in range(n_l):
        if Y.iloc[i, j]!=0j:
            m.addConstr(v_i_v_j[i,j]==l_v[i]*l_v[j])
            m.addConstr(theta_ji[i,j]==l_a[j]-l_a[i])
            m.addConstr(theta_ij[i,j]==l_a[i]-l_a[j])

            # theta[i,j]=theta_ij+theta_j-theta_i
            m.addConstr(theta[i,j]==cmath.polar(Y.iloc[i, j])[1]+theta_ji[i,j])

            # cos(theta_ij+theta_j-theta_i)
            m.addGenConstrCos(theta[i,j], cos_theta[i,j])

            # sin(theta_ij+theta_j-theta_i)
            m.addGenConstrSin(theta[i,j], sin_theta[i,j])

            # obj3=....2*cos(theta_j-theta_i)         
            m.addGenConstrCos(theta_ji[i,j], cos_theta_ji[i,j])

    # ************************************************
    # active power flow equation
    # Y_ij*v_j*v_i*cos(theta_ij+theta_j-theta_i)
    m.addConstr(inject_active-load.iloc[i, 0]==quicksum(
        [cmath.polar(Y.iloc[i, j])[0]*v_i_v_j[i, j]*cos_theta[
            i,j] for j in range(n_l) if Y.iloc[i, j]!=0j]
        ))
    # ************************************************
    # reactive power flow equation
    # -Y_ij*v_j*v_i*sin(theta_ij+theta_j-theta_i)
    m.addConstr(inject_reactive-load.iloc[i, 1]==-quicksum(
        [cmath.polar(Y.iloc[i, j])[0]*v_i_v_j[i, j]*sin_theta[
            i,j] for j in range(n_l) if Y.iloc[i, j]!=0j]
        ))

params = [0.0045, 87, 25, 92, 92]


# *********** Set objective
m.update()

obj1 = params[0]*g_p**2+params[1]*g_p+params[2]
obj2 = params[3]*s[0]

obj3 = params[4]*1/2*quicksum(
    [cmath.polar(Y.iloc[i, j])[0]*np.cos(
        cmath.polar(Y.iloc[i, j])[1])*(l_v[i]**2+l_v[j]**2-v_i_v_j[
            i, j]*2*cos_theta_ji[i,j]) for i in range(
                n_l) for j in range(n_l) if Y.iloc[i, j]!=0j])

m.setObjective(obj1+obj2+obj3, GRB.MINIMIZE) #

m.optimize()

# Set some parameters, like TIME_LIMIT, to determine 
# the termination before obtaining the optimal solution. 
# if m.status == GRB.Status.OPTIMAL:
# else:
#     print("******************************************************")
#     print("There is no feasible solution")
  
for i in m.getVars()[:2+2+n_l*2]:
    print(i.varName, i.X)
solu = [var.X for var in m.getVars()]


print('running time:', time.time()-start_time)

# solution analysis
print("active generation", solu[0]+solu[2])
print("reactive generation", solu[1]+solu[3])
print("active load", load.sum()[0])
print("reactive load", load.sum()[1])
