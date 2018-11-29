
import numpy as np
from scipy import integrate as igt
import time

import math
# def cks(t,Y):
#     Nsp = 13
#     T = Y[Nsp]
#     d = np.empty([Nsp+1, 1])
#     for i in range(Nsp):
#         k[i] = A[i]*T**b[i]*math.exp(-E[i]/R/T)
#         s[i] = h[i]*Y[i]*odot[i]
#         d[i] = Y[i]/rho*odot[i]
#
#     d[Nsp] = -1 / rho / cp * fsum(s)
#
#     return d
#
# Y0 = np.empty([15, 1])
# Y0[0] = 9.07*10**-12
# Y0[1] = 0.02789
# Y0[2] = 1.92*10**-11
# Y0[3] = 1.25*10**-10
# Y0[4] = 0.005443
# Y0[5] = 0.2212
# Y0[6] = 7.78*10**-6
# Y0[7] = 2.65*10**-4
# Y0[8] = 0
# Y0[9] = 0
# Y0[10] = 0
# Y0[11] = 0
# Y0[12] = 0
# Y0[13] = 0.7451
# Y0[14] = 850.48



#################################################################

def simple_cks(t,y):

    k1 = 2 * 10 ** 3
    k2 = 1 * 10 ** -3
    k3 = 10
    Nsp = 4
    dY = np.empty([Nsp, 1])
    dY[0] = -k1 * y.item(0) * y.item(1) + k2 * y.item(2)
    dY[1] = -k1 * y.item(0) * y.item(1) + (k2 + k3) * y.item(2)
    dY[2] = k1 * y.item(0) * y.item(1) - (k2 + k3) * y.item(2)
    dY[3] = k3 * y.item(2)

    return dY.reshape(Nsp)


def runsolver(sol,f,Solver):
    print(Solver)
    print('Solving...')

    while sol.t < sol.t_bound:
        rt0 = time.time()
        sol.step()
        rt1 = time.time()

        s = str(sol.step_size) + " " + str(rt1-rt0) + " " + str(sol.t) + " " + str(sol.y) + " " + '\n'
        f.write(s)
        print(sol.t)
    return print('Finished!')

t0 = 0
tf = 3000
Y0 = (1, 5*10**-5, 0, 0)


Solver = 'RK4'
sol = igt.RK45(simple_cks, t0, Y0, tf)  # initialize the solution
print(sol.y)
fileName = 'simple_RK4_data.txt'
f = open(fileName,'w')
# runsolver(sol,f,Solver)
f.close()

#######################################################

Solver = 'BDF'
sol = igt.BDF(simple_cks, t0, Y0, tf)
fileName = 'simple_BDF_data.txt'
f = open(fileName,'w')
# runsolver(sol,f,Solver)
f.close()