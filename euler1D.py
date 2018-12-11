#Programmer: Anthony Walker
"""This is an implementation of the sod shock tube problem."""
#READ ME:
#This is an implementation of the sod shock tube problem
#In order to execute this implementation call the euler1D function
#Changes to the time step, domain, and other properties can be made in
#problem setup section below.

#Import Statements
import fluid_domain as fd
import numpy as np
import time as tP
import os
import math
import supportingFiles as sf
import shocktubecalc as stc

#Global variables - Allocation
gamma = 1.4
dtdx = 0
numsaves = 0
t = 0
dx = 0
saveFactor = 10 #Controls number of times to save
mG = 0
# Functions for sod shock problem
def euler1D(domain,time, g = 1.4, directory = None,**kwargs):
    """Use this method to execute euler1D."""
    #Preliminary Steps
    #Updates Needed Global Variables
    globalVariableHandler(gamma,domain.dx,time)
    #Creation of document information/heading
    eulerInfo = documentInfoGeneration()
    #Handling directory & initial domain save
    dirStr = directoryHandler(directory)
    #Determining points to save
    saves = determineSavePts()
    kwargs["saves"] = saves
    #Other Variables used for data storage
    #Temporary conversion of domain to np.array
    data = np.zeros((dims[0],len(domain.primaryDomain[0][0][:])))
    for x in range(len(data)):
        data[x,:] = domain.primaryDomain[x][0][:]
    #Calculation Begins Here
    print("Beginning Calculations...")
    RK2(fv5p,data,**kwargs)
    print("Calculation Complete...")

def RK2(fcn,data,**kwargs):
    """Use this method to solve a function with RK2 in time."""
    tCurr = 0
    sI = 1
    count = 0
    multipler = 1
    tEnd = t[0]
    dt = t[1]

    if 'save' in kwargs:
        save = True
    else:
        save = False
    #Making Q from node Data
    Q = makeQ(data[:,1:])
    #Creating pressure vector
    P = data[:,0]
    while tCurr <= tEnd:
        QS = Q+dtdx*0.5*fcn(Q,P)
        P = eqnStateQ(QS[:,0],QS[:,1],QS[:,2])
        Q = Q+dtdx*fcn(QS,P)
        if "save" in kwargs:
            print(tCurr)
            sI,count,multipler,eulerInfo,eStr = saveDataUpdater(sI,count,multipler,tCurr,eulerInfo)
            kwargs["data"] = makeND(Q)
            print(kwargs)
            input()
            dataSave(sI,count,multipler,eulerInfo,eStr,**kwargs,)
        tCurr+=dt
    return makeND(Q)

#Numerical Solution
def fv5p(Q,P):
    """Use this 5 point finite volume method to solve the euler equations."""
    #Creating flux array
    Flux = np.zeros((len(Q),len(Q[1,:])))
    #Setting up loop to span entire domain except end points
    for x in range(2,len(Q)-2):
        #Q on the left side of the minus 1/2 interface
        QLm = limiter(Q[x-1],Q[x],Q[x-1],(P[x]-P[x-1]),(P[x-1]-P[x-2]))
        #Q on the right side of the minus 1/2 interface
        QRm = limiter(Q[x],Q[x-1],Q[x],(P[x]-P[x-1]),(P[x+1]-P[x]))
        #Q on the left side of the plus 1/2 interface
        QLp = limiter(Q[x],Q[x+1],Q[x],(P[x+1]-P[x]),(P[x]-P[x-1]))
        #Q on the right side of the plus 1/2 interface
        QRp = limiter(Q[x+1],Q[x],Q[x+1],(P[x+1]-P[x]),(P[x+2]-P[x+1]))
        #Getting Flux
        Flux[x,:] += makeFlux(QLm,QRm)
        Flux[x,:] -= makeFlux(QLp,QRp)
        Flux[x,:] += spectral(QLm,QRm)
        Flux[x,:] -= spectral(QLp,QRp)
        Flux[x,:] = Flux[x,:]*0.5
    return Flux

#Step 1, make Q from node data
def makeQ(Nargs):
    """Use this method to make Q."""
    dim1 = len(Nargs)
    dim2 = len(Nargs[1,:])
    Q = np.zeros((dim1,dim2))
    k = 0
    for i in Nargs:
        Q[k,:] = np.array([i[0],i[1]*i[0],i[0]*i[2]])
        k +=1
    return Q

#Step 2, get reconstructed values of Q at the interfaces
def limiter(Q1,Q2,Q3,num,denom):
    """Use this method to apply the minmod limiter."""
    if (num > 0 and denom > 0) or (num < 0 and denom < 0):
        Q_r = Q1+min(num/denom,1)/2*(Q2-Q3)
    else:
        Q_r = Q1
    return Q_r

#Step 3, make the flux
def makeFlux(QL,QR):
    """Use this method to make Q."""
    uL = QL[1]/QL[0]
    uR = QR[1]/QR[0]
    PL = eqnStateQ(QL[0],QL[1],QL[2])
    FL = np.array([QL[1],QL[1]*uL+PL,(QL[2]+PL)*uL])
    PR = eqnStateQ(QR[0],QR[1],QR[2])
    FR = np.array([QR[1],QR[1]*uR+PR,(QR[2]+PR)*uR])
    return FL+FR

#Step 4, spectral method
def spectral(QL,QR):
    """Use this function to apply the Roe average."""
    Qsp = np.zeros((len(QL)))
    rootrhoL = np.sqrt(QL[0])
    rootrhoR = np.sqrt(QR[0])
    uL = QL[1]/QL[0]
    uR = QR[1]/QR[0]
    eL = QL[2]/QL[0]
    eR = QR[2]/QR[0]
    denom = 1/(rootrhoL+rootrhoR)
    Qsp[0] = (rootrhoL*rootrhoR)
    Qsp[1] = (rootrhoL*uL+rootrhoR*uR)*denom
    Qsp[2] = (rootrhoL*eL+rootrhoR*eR)*denom
    pSP = eqnStateQ(Qsp[0],Qsp[0]*Qsp[1],Qsp[0]*Qsp[2])
    rSP = np.sqrt(gamma*pSP/Qsp[0])+abs(Qsp[1])
    Q_rs = rSP*(QL-QR)
    return Q_rs

#The equation of state using Q variables
def eqnStateQ(r,rU,rE):
    """Use this method to solve for pressure."""
    P = mG*(rE-rU*rU/(r*2))
    return P

#Create node data
def makeND(Qargs):
    """Use this function to make node data."""
    nodeData = np.zeros((len(Qargs),len(Qargs[1,:])+1))
    P = eqnStateQ(Qargs[:,0],Qargs[:,1],Qargs[:,2])
    i = 0
    for x in Qargs:
        nodeData[i,:] = np.array([P[i],x[0],x[1]/x[0],x[2]/x[0]])
        i+=1
    return nodeData

#Analytical Solution Validation
def PyPiSol(dirStr,eStr,eulerInfo,tCurr):
    """Use this PyPi, solution to validate."""
    ss = stc.solve(t = tCurr,**{'npts':dims[0]})
    ss = ss[2]
    tDomain = tuple()
    sodStr = dirStr+"/pypiSol1"+eStr+".txt"
    for x in range(len(ss['rho'])):
        e = sf.eqnState(ss['p'][x],ss['rho'][x],ss['u'][x],gamma)
        tDomain+=(np.array([ss['p'][x],ss['rho'][x],ss['u'][x],e]),)
    sf.SolFile(tDomain,sodStr,eulerInfo)

#Analytical Solution
def analyticalSol(dirStr,eStr,eulerInfo,tCurr):
    """Use this analytical solution to compare to the numerical solution."""
    sodStr = dirStr+"/ssSol1"+eStr+".txt"
    sf.sodShock(sodStr,eulerInfo,tCurr,npts = dims[0])

#Function to save appropriate data
def dataSave(**kwargs):
    if saves[sI] == tCurr:
        if 'numerical' in kwargs:
            eulerSaveStr = dirStr+"/e1"+eStr+".txt"
            sf.SolFile(kwargs['data'],eulerSaveStr,eulerInfo)
        if 'analytical' in kwargs:
            analyticalSol(dirStr,eStr,eulerInfo,tCurr)
        if 'pypi' in kwargs:
            PyPiSol(dirStr,eStr,eulerInfo,tCurr)

#Function to update saving data
def saveDataUpdater(sI,count,multipler,tCurr,eulerInfo):
    """Use this method to update strings and data for saving."""
    print("Current time:"+s1)
    s1 = "%.4f" % tCurr
    #This is for sorting purposes/data analysis
    if(96+sI)<=122:
        eStr = chr(96+sI)
    else:
        if count == 10:
            count = 0
            multipler += 1
        eStr = chr(122)*multipler+str(count)
        count+=1
    eulerInfo[0] = ("Euler test case, time = "+s1)
    sI+=1
    return (sI,count,multipler,eulerInfo,eStr,)

#Determine points to save data
def determineSavePts():
    """Use this method to determine which time indices to save."""
    timeSteps = time[0]/time[1]
    if timeSteps%10>0:
        timeSteps = int(np.ceil(timeSteps))
    else:
        timeSteps = int(timeSteps)
    index = timeSteps/numSaves
    saveInd = tuple()
    for i in range(timeSteps+1):
        if(i%index == 0):
            saveInd+=(i,)
    return saveInd

#Creates directory if necessary
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

#Handles Global Variables
def globalVariableHandler(g,dX,time):
    """Use this method to update global variables before calculation."""
    global gamma
    gamma = g #Allows changings of gamma if =desired
    global mG
    mG = gamma -1
    #Stepsize
    global dx
    dx = dX
    #Determining time step
    global dtdx
    dtdx = time[1]/dx
    #Number of times to saves
    global numSaves
    numSaves = time[0]*saveFactor
    # Time steps
    global t
    t = time
    return

#Creates document header
def documentInfoGeneration():
    """Use this to generate documentation heading information"""
    """Use this to generate documentation heading information"""
    eulerInfo = list()
    eulerInfo.append("Euler test case, time = 0")
    eulerInfo.append("dx = "+str(dx))
    eulerInfo.append("dt = "+str(t[1]))
    eulerInfo.append("nSteps = "+str(t[0]/t[1]))
    eulerInfo.append("[Pressure]   [Density]    [Velocity]    [Internal Energy]")
    return eulerInfo

if __name__ == "__main__":
    #Problem set up - Sod Shock Problem
    #Initial and Boundary conditions
    leftBC = (1.0,1.0,0,2.5)
    rightBC = (0.1,0.125,0,2)
    #Domain Creation and Initialization
    sf.generateNodeFile("textFiles/euler1D.txt", range(0,251), range(0,1))
    domain = fd.domain("textFiles/euler1D.txt")
    dims = domain.getDomainDims()
    domain.setNodeVals(rightBC,range(int(dims[0]/2),dims[0]),range(dims[1]))
    domain.setNodeVals(leftBC,range(0,int(dims[0]/2)),range(dims[1]))
    time = (0.1,0.0001)
    kwgs = {'save',True,'numerical',True,'analytical',True,'pypi',True}
    euler1D(domain,time,kwargs = kwgs)
