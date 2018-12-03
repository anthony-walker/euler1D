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
import math
import sodShock as sd
import nodeFileGenerator as nfg
#Global variables - Allocation
gamma = 1.4
dtdx = 0
numsaves = 0
timeSteps = 0
dx = 0

# Functions for sod shock problem
def euler1D(domain,time, g = 1.4, directory = None,sA = False):
    """Use this method to execute euler1D."""
    #Preliminary Steps
    #Updates Needed Global Variables
    globalVariableHandler(gamma,domain.dx,time)
    #Creation of document information/heading
    eulerInfo = documentInfoGeneration()
    #Handling directory & initial domain save
    dirStr = directoryHandler(directory)
    domain.domainToFile(dirStr+"/e1a.txt",eulerInfo)
    #Determining points to save
    saves = determineSavePts()
    #Other Variables used for data storage
    tCurr = 0
    sI = 1
    count = 0
    multipler = 1
    #Calculation Begins Here
    print("Beginning Calculations...")
    for k in range(timeSteps+1):
        #RK2 in Time Euler Equation Calculation
        qHalfStep = euler(0)
        for i in range(len(qHalfStep)):
            domain[i+1,0][:] = qHalfStep[i]
        qStep = euler(1)
        for i in range(len(qHalfStep)):
            domain[i+1,0][:] = qHalfStep[i]
        #Saving Data
        if saves[sI] == k:
            s1 = "%.4f" % tCurr
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
            if(sA):
                sodStr = dirStr+"/aSol1"+eStr+".txt"
                sd.sodShock(sodStr,eulerInfo,tCurr,numPts = dims[0])
        tCurr+=time[1]
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
    for i in range(1,dims[0]-1):
        #Getting points from domain
        P = domain[i,j][:]
        W = boundaryHandler(domain,i,j,-1)
        WW = boundaryHandler(domain,i,j,-2)
        E = boundaryHandler(domain,i,j,1)
        EE = boundaryHandler(domain,i,j,2)
        #nodes
        ns = (WW,W,P,E,EE,)
        #makeQ
        qN = makeQ(ns)
        #Getting reconstructed Q values at interface
        qI = mmdlim(qN,0,i)
        #finding flux
        fU = flux(qI[0],qI[1])
        fD = flux(qI[2],qI[3])
        flx += unmakeQ((qN[2][1:] + dtdx*sc*(fU-fD),))
    return flx

def flux(qL,qR):
    """Use this method to determine the flux."""
    #Preallocation
    f = np.zeros(len(qL))
    qSP = np.zeros(len(qL))
    half = 0.5
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
    qSP[2] = (rootrhoL*eL+rootrhoR*eR)/(rootrhoL+rootrhoR)
    pSP = eqnState(qSP[0],qSP[1],qSP[2])
    rSP = np.sqrt(gamma*pSP/qSP[0])+abs(qSP[1])

    #flux
    pL = eqnState(qL[0],uL,eL)
    pR = eqnState(qR[0],uR,eR)
    f[0] = half*(qL[1]+qR[1]+rSP*(qL[0]-qR[0]))
    f[1] = half*(qL[0]*uL**2+pL+qR[0]*uR**2+pR+rSP*(qL[1]-qR[1]))
    f[2] = half*(uL*(qL[2]+pL)+uR*(qR[2]+pR)+rSP*(qL[2]-qR[2]))
    return f
def eqnState(rho,u,e):
    """Use this method to solve for pressure."""
    P = (gamma-1)*(e-rho*u*u/2)
    return P
def mmdlim(P,limInd,k):
    """Use this method to get Q with a flux limitor."""
    #Instance tuples
    d = tuple()
    Q = tuple()
    C = tuple()
    x = 1 #Handles exponent for pressure ratio
    c = 0 #Step counter
    a = 0 #Difference counter
    #Loop to get diffs and values of Q for minmod 0, 1, 1, 2 as example
    for i in range(len(P)-1): #Pressure differences for all 5 points
        d += (P[i+1][limInd]-P[i][limInd],)
        if(x == 1):
            x = x*-1
            C+=(c,)
            c = c+1
        else:
            x = x*-1
            C+=(c,)
    P = P[1:len(P)-1] #Adjusting P to be only W, P, E
    for c in C:
        if (d[c] < 0 and d[c-1] < 0) or (d[c] > 0 and d[c-1] > 0): #Checking Differences
            Pr = d[c]/d[c-1]
            Q += (P[c][1:]+min(Pr**x,1)/2*x*(P[a+1][1:]-P[a][1:]),)
        else:
            Q += (P[c][1:],)
        if x < 0:
            a += 1
        x = x*-1
    return Q

def unmakeQ(Q):
    """Use this method to obtain Q values from nodeValues."""
    nV = tuple()
    for x in Q:
        tA = np.array([x[0],x[1]/x[0],x[2]/x[0]])
        P = eqnState(tA[0],tA[1],tA[2])
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
    return saveInd

def directoryHandler(directory = None):
    """Use this method to check if a directory exists or make one."""
    if(directory is None):
        directory = input("Enter a directory for values to be stored inside results. ")
        dirStr = "results/"+directory
    else:
        dirStr = directory
    if not os.path.isdir(dirStr):
        print("Creating directory...")
        os.mkdir(dirStr)
    return dirStr

def globalVariableHandler(g,dX,time):
    """Use this method to update global variables before calculation."""
    global gamma
    gamma = g #Allows changings of gamma if desired
    #Stepsize
    global dx
    dx = dX
    #Determining time step
    global dtdx
    dtdx = time[1]/dx
    #Number of times to saves
    global numSaves
    numSaves = time[0]*50
    # Time steps
    global timeSteps
    timeSteps = time[0]/time[1]
    if timeSteps%10>0:
        timeSteps = int(np.ceil(timeSteps))
    else:
        timeSteps = int(timeSteps)
    return

def documentInfoGeneration():
    """Use this to generate documentation heading information"""
    eulerInfo = list()
    eulerInfo.append("Euler test case, time = 0")
    eulerInfo.append("dx = "+str(dx))
    eulerInfo.append("dt = "+str(time[1]))
    eulerInfo.append("nSteps = "+str(timeSteps))
    eulerInfo.append("[Pressure]   [Density]    [Velocity]    [Internal Energy]")
    return eulerInfo

if __name__ == "__main__":
    #Problem set up - Sod Shock Problem
    #Initial and Boundary conditions
    leftBC = (1.0,1.0,0,2.5)
    rightBC = (0.1,0.125,0,0.25)
    #Domain Creation and Initialization
    nfg.generateNodeFile("euler1D.txt", range(0,1001), range(0,1))
    domain = fd.domain("euler1D.txt")
    dims = domain.getDomainDims()
    domain.setNodeVals(rightBC,range(int(dims[0]/2),dims[0]),range(dims[1]))
    domain.setNodeVals(leftBC,range(0,int(dims[0]/2)),range(dims[1]))
    time = (1,0.000001)
    euler1D(domain,time,sA = True)
