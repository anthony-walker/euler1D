#Programmer: Anthony Walker
#This is a five point finite volume method used to solve the euler equations

import numpy as np
import supportingFiles as sf

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
        P = eqnState(QS[:,0],QS[:,1]/QS[:,0],QS[:,2]/QS[:,0])
        Q = Q+dt/dx*fcn(QS,P)
        tCurr+=dt
        print(makeND(Q))
        input()
    return makeND(Q)

def fv5p(Q,P):
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
        Flux[x-2,:] += makeFlux(QLm,QRm)
        Flux[x-2,:] -= makeFlux(QLp,QRp)
        Flux[x-2,:] += spectral(QLm,QRm)
        Flux[x-2,:] -= spectral(QLp,QRp)
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
    PL = eqnState(QL[0],QL[1]/QL[0],QL[2]/QL[0])
    FL = np.array([QL[1],QL[1]*QL[0]+PL,(QL[2]+PL)*QL[1]/QL[0]])
    PR = eqnState(QR[0],QR[1]/QR[0],QR[2]/QR[0])
    FR = np.array([QR[1],QR[1]*QR[0]+PR,(QR[2]+PR)*QR[1]/QR[0]])
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
    pSP = eqnState(Qsp[0],Qsp[1],Qsp[2])
    rSP = np.sqrt(gamma*pSP/Qsp[0])+abs(Qsp[1])
    Q_rs = rSP*(QL-QR)
    return Q_rs

def eqnState(rho,u,e):
    """Use this method to solve for pressure."""
    P = mG*(e-rho*u*u/2)
    return P

def makeND(Qargs):
    """Use this function to make node data."""
    nodeData = np.zeros((len(Qargs),len(Qargs[1,:])+1))
    P = eqnState(Qargs[:,0],Qargs[:,1]/Qargs[:,0],Qargs[:,2]/Qargs[:,0])
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
            currArray[x,:] = np.array([0.1,0.125,0,0.25])
    print("Begin")
    eulerInfo = documentInfoGeneration()
    currArray = RK2(fv5p,currArray,dt,dx,t)
    sf.analytSolFile(currArray,"/home/walkanth/euler1D/results/test/newNumSol.txt",eulerInfo)
    print("End")
