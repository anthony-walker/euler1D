"""
Created on Tues Nov 20 14:54:33 2018

@author: Paige
"""


import numpy as np
from scipy import integrate as igt
import time
import math


def oregonator(t,y):
    s = 77.27
    q = 8.375*10**-6
    w = 0.161

    dY = np.empty([3, 1])
    dY[0] = s*(y.item(0) - y.item(0)*y.item(1) + y.item(1) - q*y.item(0)*y.item(0))
    dY[1] = (y.item(2) - y.item(1) - y.item(0)*y.item(1)) / s
    dY[2] = w*(y.item(0)-y.item(2))

    return dY.reshape(3)

def reinitialsolver(sol,f,Solver,globaltimestep):
    print(Solver)
    print('Solving...')

    n = math.floor(sol.t_bound/globaltimestep)
    r = sol.t_bound/globaltimestep - n

    for i in range(n):

        if i is n:
            t_b = (n+r) * globaltimestep
        else:
            t_b = (i+1) * globaltimestep
        # print(t_b)
        t0 = sol.t
        Y0 = sol.y

        if Solver is 'RK4':
            sol = igt.RK45(oregonator, t0, Y0, t_b)
        if Solver is 'BDF':
            sol = igt.BDF(oregonator, t0, Y0, t_b)
        grt0 = time.time()
        while sol.t < sol.t_bound:

            rt0 = time.time()
            sol.step()
            rt1 = time.time()

            s = 'NaN '+str(sol.step_size) + " " + str(rt1-rt0) + " " + str(sol.t) + " " + str(sol.y) + " " + '\n'
            f.write(s)
            # print(sol.t)
        grt1 = time.time()
        f.write(str(grt1-grt0)+' NaN NaN NaN NaN \n')

    return print('Finished!')


def runsolver(sol,f,Solver):
    print(Solver)
    print('Solving...')

    while sol.t < sol.t_bound:
        rt0 = time.time()
        sol.step()
        rt1 = time.time()

        s = str(sol.step_size) + " " + str(rt1-rt0) + " " + str(sol.t) + " " + str(sol.y) + " " + '\n'
        f.write(s)
        # print(sol.t)
    return print('Finished!')

t0 = 0
tf = 375
Y0 = (1, 1, 2)
p = 1
gts = 10**-p


Solver = 'RK4'
sol = igt.RK45(oregonator, t0, Y0, tf)  # initialize the solution
sol.atol = 10**-10
sol.rtol = 10**-6
fileName = 'oregonator_RK4_data_En'+ str(p) +'.txt'
f = open(fileName,'w')
reinitialsolver(sol,f,Solver,gts)
f.close()

Solver = 'BDF'
sol = igt.BDF(oregonator, t0, Y0, tf)  # initialize the solution
sol.atol = 10**-10
sol.rtol = 10**-6
fileName = 'oregonator_BDF_data_En'+ str(p) +'.txt'
f = open(fileName,'w')
reinitialsolver(sol,f,Solver,gts)
f.close()