# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 18:56:33 2018

@author: Paige
"""

import numpy as np
from scipy import integrate as igt




def vanderPol(t, y):
    # y = np.empty([2,1])
    u = 1000.
    dy = [y[1], -y[0]+u*y[1]*(1-y[0]**2)]
    return dy


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
y0 = (2, 0)
# y0 = (1.73375132e+00, -8.62850136e-04)

t_bound = 400
# dt_min = 0.00035
dydt_crit = 0.0007
y = list(y0)
t = list([t0])

sol = igt.RK45(vanderPol, t0, y0, t_bound)  # initialize the solution
# sol.step()
# t0 = sol.t
# y0 = sol.y
k = bool(1)
print(type(k))
print(k)
while sol.t < t_bound:
    if abs(round(sol.y[1], 7)) < dydt_crit and k is True:
        print('Solving Non-Stiff Portion...')
        print(sol.step_size, ',  ', sol.y, ',  ', sol.t)
        t0 = sol.t
        y0 = sol.y
        sol = igt.RK45(vanderPol, t0, y0, t_bound)
        k = bool(0)

    if abs(round(sol.y[1], 7)) >= dydt_crit*0.9 and k is False:
        print('Solving Stiff Portion...')
        print(sol.step_size, ',  ', sol.y, ',  ', sol.t)
        t0 = sol.t
        y0 = sol.y
        sol = igt.BDF(vanderPol, t0, y0, t_bound)
        k = bool(1)

    sol.step()
    y.append(sol.y)
    dt = sol.step_size
    t.append(sol.t)
    # print(sol.step_size, ',  ', sol.y, ',  ', sol.t)


# Data for plotting


import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(t, y)

ax.set(xlabel='time (s)', ylabel='voltage (mV)',
       title='About as simple as it gets, folks')
ax.grid()
plt.show()














