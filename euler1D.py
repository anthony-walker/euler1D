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
import time as tm
import os
import math
import supportingFiles as sf
import shocktubecalc as stc
import datetime as dt

#Global Variables
startTime = tm.time() #run time start
print("Starting time: ", str(dt.datetime.now()))
#Global variables - Allocation
gamma = 1.4 #heat capacity ratio
dtdx = 0 #dt/dx
numsaves = 0 #Number of saves
t = 0 #End time and time step size
dx = 0 #Step SizeStarting
saveFactor = 10 #Controls number of times to save
mG = 0 #Gamma -1

# Functions for sod shock problem
def euler1D(domain,time, g = 1.4, directory = None,**kwargs):
    """Use this method to execute euler1D."""
    #Updates Needed Global Variables
    globalVariableHandler(gamma,domain.dx,time)
    #Creation of document information/heading
    eulerInfo = documentInfoGeneration()
    #Determining points to save
    saves = determineSavePts()
    #Handling directory & initial domain save
    dirStr = directoryHandler(directory)

    kwargs['saves'] = saves
    kwargs['eulerInfo'] = eulerInfo
    kwargs['dirStr'] = dirStr
    #Other Variables used for data storage
    #Temporary conversion of domain to np.array
    data = np.zeros((dims[0],len(domain.primaryDomain[0][0][:])))
    for x in range(len(data)):
        data[x,:] = domain.primaryDomain[x][0][:]
    #Calculation Begins Here
    print("Beginning Calculations...")
    RK2(fv5p,data,**kwargs)
    print("Calculation Complete...")
    print("Runtime: ",(tm.time()-startTime)/60," minutes.")
    print("Ending time: ", str(dt.datetime.now()))
    askToPlot(dirStr+"/",eulerInfo)

### Functions Crucial To Calculation ###
def RK2(fcn,data,**kwargs):
    """Use this method to solve a function with RK2 in time."""
    #Preliminary Steps
    tCurr = 0
    sI = 1
    count = 0
    multipler = 1
    tEnd = t[0]
    dt = t[1]
    step = 0
    #Making Q from node Data
    Q = makeQ(data[:,1:])
    #Creating pressure vector
    P = data[:,0]
    while tCurr <= tEnd:
        QS = Q+dtdx*0.5*fcn(Q,P)
        P = eqnStateQ(QS[:,0],QS[:,1],QS[:,2])
        Q = Q+dtdx*fcn(QS,P)
        tCurr+=dt
        step+=1
        #Saving data if desired
        if sI < len(kwargs['saves']): #Extra loop as a precaution
            if kwargs["save"] and kwargs['saves'][sI]==step:
                sI,count,multipler,kwargs['eulerInfo'],kwargs["eStr"] = saveDataUpdater(sI,count,multipler,tCurr,kwargs['eulerInfo'])
                kwargs["data"] = makeND(Q)
                dataSave(tCurr,**kwargs)
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
        nodeData[i,:] = np.array([P[i],x[0],x[1]/x[0],x[2]])
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

### Functions Not Crucial To Calculation ###
#Function to save appropriate data
def dataSave(tCurr,**kwargs):
    """Use this method to save all data."""
    dirStr = kwargs["dirStr"]
    eulerInfo = kwargs["eulerInfo"]
    eStr = kwargs['eStr']
    if kwargs["numerical"]:
        eulerSaveStr = dirStr+"/e1"+eStr+".txt"
        sf.SolFile(kwargs["data"],eulerSaveStr,eulerInfo)
    if kwargs["analytical"]:
        analyticalSol(dirStr,eStr,eulerInfo,tCurr)
    if kwargs["pypi"]:
            PyPiSol(dirStr,eStr,eulerInfo,tCurr)

#Function to update saving data
def saveDataUpdater(sI,count,multipler,tCurr,eulerInfo):
    """Use this method to update strings and data for saving."""
    s1 = "%.4f" % tCurr
    print("Saving Current time:"+s1,"...")
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
    print("Determining times to save...")
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
    print("Handling Global Variables...")
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
    print("Generating Document Information...")
    eulerInfo = list()
    eulerInfo.append("Euler test case, time = 0")
    eulerInfo.append("dx = "+str(dx))
    eulerInfo.append("dt = "+str(t[1]))
    eulerInfo.append("nSteps = "+str(t[0]/t[1]))
    eulerInfo.append("[Pressure]   [Density]    [Velocity]    [Internal Energy]")
    return eulerInfo

def askToPlot(dirStr,eulerInfo):
    """Use this function to ask if the user wants to plot data."""
    pB = input("Would you like to plot the data (Y/N)? ")
    if pB is 'Y':
        labels = ["Pressure","Density","Velocity","Internal Energy"]
        lSkip = len(eulerInfo)+1
        sf.dataPlot(labels,lSkip,dirStr)

if __name__ == "__main__":
    #Problem set up - Sod Shock Problem
    #Initial and Boundary conditions
    leftBC = (1.0,1.0,0,2.5)
    rightBC = (0.1,0.125,0,2)
    #Domain Creation and Initialization
    sf.generateNodeFile("textFiles/euler1D.txt", range(0,101), range(0,1))
    domain = fd.domain("textFiles/euler1D.txt")
    dims = domain.getDomainDims()
    domain.setNodeVals(rightBC,range(int(dims[0]/2),dims[0]),range(dims[1]))
    domain.setNodeVals(leftBC,range(0,int(dims[0]/2)),range(dims[1]))
    time = (0.3,0.001)
    print("time end: ",time[0]," dt: ",time[1])
    print("dx: ", domain.dx)
    kwgs = {"save":True,"numerical":True,"analytical":True,"pypi":False}
    euler1D(domain,time,**kwgs)
