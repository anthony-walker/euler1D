#Testing script
import shocktubecalc as stc
import numpy as np
import datetime as dt

ss = stc.solve(t = 0.2,**{'npts':251})
ss = ss[2]
gamma = 1.4
tDomain = tuple()
fileName = 'results/testingA.txt'
dx = 1/250

"""Use this to generate documentation heading information"""
eulerInfo = list()
eulerInfo.append("Euler test case, time = 0")
eulerInfo.append("dx = "+str(0))
eulerInfo.append("dt = "+str(0))
eulerInfo.append("nSteps = "+str(0))
eulerInfo.append("[Pressure]   [Density]    [Velocity]    [Internal Energy]")
extraInfo = eulerInfo

for x in range(len(ss['rho'])):
    e = 1#eqnState(ss['p'][x],ss['rho'][x],ss['u'][x],gamma)
    tDomain+=(np.array([ss['p'][x],ss['rho'][x],ss['u'][x],e]),)
if(fileName is None):
    fileName = input("Enter file name and extension. ")
data = tDomain
with open(fileName, 'w') as f:

    f.write("Date & Time: ")
    f.write(str(dt.datetime.now()))
    f.write("\n")

    if extraInfo is not None:
        for e in extraInfo:
            f.write(e)
            f.write("\n")
    for x in data:
        for y in x:
            f.write(" %0.8f " % y)
        f.write("\n")
    f.closed
