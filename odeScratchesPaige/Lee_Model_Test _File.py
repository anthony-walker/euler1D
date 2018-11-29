
"""
Created on Mon Nov  12 15:53:33 2018

@author: Paige
"""


import numpy as np
from scipy import integrate as igt
import time



def leeModel(t, y):
    e = 10. ** -2
    Nsp = 4
    dy = np.empty([Nsp, 1])
    s = 0
    print(y)
    for j in range(Nsp - 2):
        s = + y.item(j + 1) / (1 + y.item(j + 1))

    for i in range(Nsp - 2):
        dy[i] = [1 / e ** (Nsp - i) * (-y.item(i) + y.item(i + 1) / (1 + y.item(i + 1)) - s)]

    dy[Nsp - 1] = -y.item(Nsp - 1)

    return dy.reshape(Nsp)

def runsolver(sol,f,Solver):
    print(Solver)
    print('Solving...')

    while sol.t < sol.t_bound:
        if Solver is 'RK4':
            rt0 = time.time()
            sol.step()
            sol.step()
            sol.step()
            sol.step()
            sol.step()
            sol.step()
            sol.step()
            sol.step()
            sol.step()
            sol.step()
            rt1 = time.time()
        else:
            rt0 = time.time()
            sol.step()
            rt1 = time.time()

        s = str(sol.step_size) + " " + str(rt1-rt0) + " " + str(sol.t) + " " + str(sol.y) + " " + '\n'
        f.write(s)
        print(sol.t)
    return print('Finished!')

t0 = 0
y0 = (1, 1, 1, 1)
t_bound = 10


Solver = 'RK4'
sol = igt.RK45(leeModel, t0, y0, t_bound)  # initialize the solution
fileName = 'leeModel_RK4_data.txt'
f = open(fileName,'w')
runsolver(sol,f,Solver)
f.close()

#######################################################

Solver = 'BDF'
sol = igt.BDF(leeModel, t0, y0, t_bound)
fileName = 'leeModel_BDF_data.txt'
f = open(fileName,'w')
runsolver(sol,f,Solver)
f.close()

