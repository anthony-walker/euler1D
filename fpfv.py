#Programmer: Anthony Walker
#This is a five point finite volume method used to solve the euler equations

import numpy as np
import supportingFiles as sf
import time
#Temporary global variables
gamma = 1.4
mG = gamma-1
def RK2(fcn,data,dt,dx,t):
    """Use this method to solve a function with RK2 in time."""
    tCurr = 0
    #Making Q from node Data
    Q = makeQ(data[:,1:])
    #Creating pressure vector
    P = data[:,0]
    while tCurr <= t:
        QS = Q+dt/dx*0.5*fcn(Q,P)
        P = eqnStateQ(QS[:,0],QS[:,1],QS[:,2])
        Q = Q+dt/dx*fcn(QS,P)
        tCurr+=dt
    return makeND(Q)

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
    return Q_rsdef RK2(fcn,domain):
    """Use this method to apply RK2 in time."""
    f = fd.domain()
    f[:] = domain #Original domain
    qHalfStep = fcn(domain)
    for i in range(len(qHalfStep)):
        qHS = makeQ((domain[i+1,0][:],))
        qHS = qHS[1:]+dtdx*qHalfStep[i]*0.5
        nVS = unmakeQ((qHS,))
        domain[i+1,0][:] = nVS
    qStep = fcn(domain)
    for i in range(len(qStep)):
        qS = makeQ((domain[i+1,0][:],))
        qS = qS[1:]+dtdx*qStep[i]
        nVS = unmakeQ((qS,))
        domain[i+1,0][:] = nVS
    return domain

def flux(qL,qR):
    """Use this method to determine the flux."""
    #Preallocation
    f = np.zeros(len(qL))
    #getting principal node values
    uL = qL[1]/qL[0]
    uR = qR[1]/qR[0]
    #flux
    pL = eqnState(qL[0],uL,qL[2]/qL[0])
    pR = eqnState(qR[0],uR,qR[2]/qR[0])
    f[0] = qL[1]+qR[1]
    f[1] = qL[1]*uL+pL+qR[1]*uR+pR
    f[2] = qL[2]*uL+pL*uL+qR[2]*uR+pR*uR
    return f

def fpfv(domain):
    """Use this to solve euler1D."""
    j=0 #This is for 2D
    flx = tuple()
    global tSum
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
        flx += (0.5*(flux(qI[0],qI[1])+spectral(qI[0],qI[1])
              -flux(qI[2],qI[3])-spectral(qI[2],qI[3])),)
    return flx

def spectral(qL,qR):
    """Use this method to calculate the spectral values"""
    #Preallocation
    qSP = np.zeros(len(qL))
    #spectral
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

def documentInfoGeneration():
    """Use this to generate documentation heading information"""
    eulerInfo = list()
    eulerInfo.append("Euler test case, time = 0")
    eulerInfo.append("dx = "+str(1/100))
    eulerInfo.append("dt = "+str(0.1))
    eulerInfo.append("nSteps = "+str(1000))
    eulerInfo.append("[Pressure]   [Density]    [Velocity]    [Internal Energy]")
    return eulerInfo

if __name__ == '__main__':
    nP = 251
    dx = 1/(nP-1)
    dt = 0.0001
    t = 0.1
    currArray = np.zeros((nP,4))
    for x in range(len(currArray)):
        if x < len(currArray)/2:
            currArray[x,:] = np.array([1,1,0,2.5])
        else:
            currArray[x,:] = np.array([0.1,0.125,0,2])
    print("Begin")
    eulerInfo = documentInfoGeneration()
    currArray = RK2(fv5p,currArray,dt,dx,t)
    # sf.analytSolFile(currArray,"/home/walkanth/euler1D/results/test/newNumSol.txt",eulerInfo)
    print("End")
