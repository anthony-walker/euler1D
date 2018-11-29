# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 18:56:33 2018

@author: Paige
"""

import numpy as np
from scipy import integrate as igt




# def vanderPol(t, y):
#     # y = np.empty([2,1])
#     u = 1000.
#     dy = [y[1], -y[0]+u*y[1]*(1-y[0]**2)]
#     return dy


def leeModel(t, y):
    print("Here")
    e = 10.**-2
    Nsp = 4
    dy = np.array([Nsp, 1])
    s = np.array(1)
    for i in range(Nsp-2):
        for j in range(i, Nsp-2):
            s[0] = sum(y[i+1]/91+y[i+1])

        dy[i] = 1/e**(Nsp-i)*(-y[i]+y[i+1]/(1+y[i+1])-s)

    dy[Nsp-1] = -y[Nsp-1]

    return dy


t0 = 0
# y0 = (2, 0)
y0 = (1.73375132e+00, -8.62850136e-04)

t_bound = 500
dt_min = 0.0004
dyDt_critical = 0.0007
y = list([y0[0]])
t = list([t0])
tol = 2*10**-5

roots = []
roots.append(abs(complex(-2.7829356341, 0)))
roots.append(abs(complex(-0.6073532183, -2.8718997282)))
roots.append(abs(complex(-0.6073532183,  2.8718997282)))
roots.append(abs(complex(-2.2194468920,  1.6872945499)))
roots.append(abs(complex(-2.2194468920, -1.6872945499)))
roots.append(abs(complex( 0.2194468920,  2.4753057481)))
roots.append(abs(complex( 0.2194468920, -2.4753057481)))

sol = igt.RK45(vanderPol, t0, y0, t_bound)  # initialize the solution
sol.step()
t0 = sol.t
y0 = sol.y
y.append(y0)
t.append(t0)
dt = sol.step_size

k = True
i = 1
n = []
Solver = []
print('Solving... ')

while sol.t < t_bound:

    if dt < dt_min and k is True:
        print('Solving Non-Stiff Portion...')
        # print(sol.step_size, ',  ', sol.y, ',  ', sol.t)
        t0 = sol.t
        y0 = sol.y
        sol = igt.RK45(vanderPol, t0, y0, t_bound)
        k = False
        i = 20
    # print(dt)
    if dt >= dt_min and k is False and i is 0:
        print('Solving Stiff Portion...')
        # print(sol.step_size, ',  ', sol.y, ',  ', sol.t)
        t0 = sol.t
        y0 = sol.y
        sol = igt.BDF(vanderPol, t0, y0, t_bound)
        k = True

    sol.step()
    y.append(sol.y[0])
    dt = sol.step_size
    t.append(sol.t)
    n.append(sol.njev)
    stepSize = sol.step_size

    if i > 0:
        i = -1

    if k:
        Solver.append(1)
    else:
        Solver.append(0)
    # print(Solver[-1])
    if sol.step_size >= dt_min:
        n.append(sol.njev)
        # print(n[-1])
        ns =+ 1

        if n[-1] is 1:
            eigValue = np.linalg.eig(sol.J)
            dt = 0.8*roots*max(abs(np.array(eigValue)))**-1

    if dt < dt_min:
        dt = sol.step_size
        s =+1


print('done')



    # print(sol.step_size, ',  ', sol.y, ',  ', sol.t)
#
# print(n[1])
#
# print('Plotting')
# y = np.array(y)
# t = np.array(t)
# # Data for plotting
#
#
# import matplotlib.pyplot as plt
#
# fig, ax = plt.subplots()
# # x = [1,2,3,4,3,2,1]
# ax.plot(t, y)
#
# ax.set(xlabel='time (s)', ylabel='voltage (mV)',
#        title='About as simple as it gets, folks')
# ax.grid()
# plt.show()














