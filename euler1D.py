#Programmer: Anthony Walker
"""This is an implementation of the sod shock tube problem."""
#READ ME:
#This is an implementation of the sod shock tube problem
#In order to execute this implementation call the euler1D function
#Changes to the time step, domain, and other properties can be made in
#problem setup section below.

import fluid_domain as fd
import numpy as np
import time
import os

#Problem set up - Sod Shock Problem
gamma = 1.4
leftBC = (1.0,1.0,0,2.5)
rightBC = (0.1,0.125,0,0.25)
domain = fd.domain("euler1D.txt")
dims = domain.getDomainDims()
domain.setNodeVals(rightBC,range(dims[0]),range(dims[1]))
domain.setNodeVals(leftBC,range(1),range(dims[1]))
dx = 1/dims[0]
dt = 0.00001
t = 2
dtdx = dt/dx
# Time steps
timeSteps = t/dt
if timeSteps%10>0:
    timeSteps = int(np.ceil(timeSteps))
else:
    timeSteps = int(timeSteps)
numSaves = t*10

# Functions for sod shock problem
def euler1D(directory = None):
    """Use this method to execute euler1D."""
    #Handling directory
    if(directory is None):
        directory = input("Enter a directory for values to be stored inside results. ")
        dirStr = "results/"+directory
    if not os.path.isdir(dirStr):
        print("Creating directory...")
        os.mkdir(dirStr)
    #Determining points to save
    saves = determineSavePts()
    eulerInfo = list()
    tCurr = 0
    sI = 1
    count = 0
    multipler = 1
    #Save Document information
    eulerInfo.append("Euler test case, time = 0")
    eulerInfo.append("dx = "+str(dx))
    eulerInfo.append("dt = "+str(dt))
    eulerInfo.append("nSteps = "+str(timeSteps))
    eulerInfo.append("[Pressure]   [Density]    [Velocity]    [Internal Energy]")
    domain.domainToFile(dirStr+"/e1a.txt",eulerInfo)

    print("Beginning Calculations...")
    for k in range(timeSteps+1):
        qHalfStep = euler(0)
        for i in range(len(qHalfStep)):
            domain[i,0][:] = qHalfStep[i]
        qStep = euler(1)
        for i in range(len(qHalfStep)):
            domain[i,0][:] = qHalfStep[i]
        if saves[sI] == k:
            s1 = "%.2f" % tCurr
            print("Current time:"+s1)
            #This is for sorting purposes/data analysis
            if(97+sI)<=122:
                eStr = chr(97+sI)
            else:
                if count == 10:
                    count = 0
                    multipler += 1
                eStr = chr(122)*multipler+str(count)
                count+=1

            eulerInfo[0] = ("Euler test case, time = "+s1)
            eulerSaveStr = dirStr+"/e1"+eStr+".txt"
            sI+=1
            domain.domainToFile(eulerSaveStr,eulerInfo)
            #domain.printNodes()
        tCurr+=dt
    print("Calculation Complete...")

def euler(step):
    """Use this to solve euler1D."""
    j=0 #This is for 2D
    flx = tuple()

    if step:
        sc = 1/2
        step = 0
    else:
        sc = 1/4
        step = 1

    for i in range(dims[0]):
        #Getting points from domain
        P = domain[i,j][:]
        W = boundaryHandler(domain,i,j,-1)
        WW = boundaryHandler(domain,i,j,-2)
        E = boundaryHandler(domain,i,j,1)
        EE = boundaryHandler(domain,i,j,2)
        #nodes
        ns = (P,W,WW,E,EE,)
        #makeQ
        qN = makeQ(ns)
        #allocate Q for interface
        qI = tuple()
        #Getting reconstructed Q values at interface
        qI += (minmod(qN[0],qN[1],qN[2],1,qN[0]-qN[1]),)
        qI += (minmod(qN[3],qN[0],qN[1],-1,qN[1]-qN[0]),)
        qI += (minmod(qN[3],qN[0],qN[1],1,qN[3]-qN[0]),)
        qI += (minmod(qN[4],qN[3],qN[0],-1,qN[0]-qN[3]),)
        #finding flux
        fU = flux(qI[0],qI[1])
        fD = flux(qI[2],qI[3])

        flx += unmakeQ((qN[0][1:] + dtdx*sc*(fU-fD),))

    return flx
def flux(qL,qR):
    """Use this method to determine the flux."""
    #Preallocation
    f = np.zeros(len(qL))
    qSP = np.zeros(len(qL))
    #getting principal node values
    rootrhoL = np.sqrt(qL[0])
    rootrhoR = np.sqrt(qR[0])
    uL = qL[1]/qL[0]
    uR = qR[1]/qR[0]
    eL = qL[2]/qL[0]
    eR = qR[2]/qR[0]
    #spectral
    qSP[0] = rootrhoL*rootrhoR
    qSP[1] = (rootrhoL*uL+rootrhoR*uR)/(rootrhoL+rootrhoR)
    qSP[2] = (rootrhoL*uL+rootrhoR*uR)/(rootrhoL+rootrhoR)
    pSP = eqnState(qSP[2],qSP[0],qSP[1])
    rSP = np.sqrt(gamma*pSP/qSP[0])+abs(qSP[1])
    #flux
    pL = eqnState(eL,qL[0],uL)
    pR = eqnState(eR,qR[0],uR)
    f[0] = 0.5*(qL[1]+qR[1]+rSP*(qL[0]-qR[0]))
    f[1] = 0.5*(qL[0]*uL**2+pL+qR[0]*uR**2+pR+rSP*(qL[1]-qR[1]))
    f[2] = 0.5*(uL*(qL[2]+pL)+uR*(qR[2]+pR)+rSP*(qL[2]-qR[2]))
    return f
def eqnState(e,rho,u):
    """Use this method to solve for pressure."""
    P = (gamma-1)*(e-rho*u**2/2)
    return P
def minmod(P1,P2,P3,expnt,diff):
    """Use this method to get Q with a flux limitor."""
    num = P1[0]-P2[0]
    den = P2[0]-P3[0]
    if(num < 0 and den < 0 or num > 0 and den > 0):
        #print("Inside minmod")
        Pr = num/den
        Q = P2[1:]+min(Pr**expnt,1)/2*diff[1:]
    else:
        #print("Inside else")
        Q = P2[1:]
    return Q
def unmakeQ(Q):
    """Use this method to obtain Q values from nodeValues."""
    nV = tuple()
    for x in Q:
        tA = np.array([x[0],x[1]/x[0],x[2]/x[0]])
        P = eqnState(tA[2],tA[0],tA[1])
        tA = np.insert(tA,0,P)
        nV += (tA,)
    return nV

def makeQ(nodes):
    """Use this method to obtain Q values from nodeValues."""
    Q = tuple()
    for x in nodes:
        Q += (np.array([x[0],x[1],x[1]*x[2],x[1]*x[3]]),)
    return Q

def boundaryHandler(pDomain,indi,indj,addition):
    """Use this method to obtain point values."""
    try:
        if (indi+addition<0):
            raise Exception("Boundary Value encountered")
        else:
            P = pDomain[indi+addition,indj][:]
    except Exception as e:
        #print(e,indi-addition,indi+addition)
        P = pDomain[indi-addition,indj][:]
    return P

def determineSavePts():
    """Use this method to determine which time indices to save."""
    index = timeSteps/numSaves
    saveInd = tuple()
    for i in range(timeSteps+1):
        if(i%index == 0):
            saveInd+=(i,)
            #print(i)
    return saveInd
if __name__ == "__main__":
    euler1D()
