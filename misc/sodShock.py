#Programmer: Anthony Walker
import numpy as np
import datetime as dt

def sodShock(fileName,fileStoreInfo,t = 0.2,lBC = (1.0,1.0,0,2.5),rBC = (0.1,0.125,0,0.25),gamma = 1.4,npts = 1000):
    """Use this method to determine the analytic solution of sod shock problem."""
    #Important coeffcients
    Gamma = (gamma-1)/(gamma+1)
    beta = (gamma-1)/(2*gamma)
    #Speed of sound for states 1 and 5
    cL = np.sqrt(gamma*lBC[0]/lBC[1])
    cR = np.sqrt(gamma*rBC[0]/rBC[1])

    #Determining states
    statesList = states(Gamma,gamma,beta,lBC,rBC,cL,cR)
    rho = statesList[0]
    u = statesList[1]
    p = statesList[2]

    #Speed of sound for state 3
    c3 = cL-(gamma-1)/2*u[2]

    #Location
    x1 = 0
    x2 = x1-cL*t
    x3 = x1 - (u[2]-c3)*t
    x4 = x1 + u[2]*t
    x5 = x1 + u[3]*t
    x = [x1,x2,x3,x4,x5]

    #Getting Data
    xI = np.linspace(-0.5,0.5,npts)
    pD = np.zeros(npts)
    rhoD = np.zeros(npts)
    uD = np.zeros(npts)
    combList = list()
    for i in range(npts):
        if(xI[i]<x1):
            rhoD[i] = rho[0]
            pD[i] = p[0]
            uD[i] = u[0]
        elif(xI[i] < x2 and xI[i] > x1):
            rhoD[i] = rho[1]
            pD[i] = p[1]
            uD[i] = u[1]
        elif(xI[i] < x3 and xI[i] > x2):
            rhoD[i] = rho[2]
            pD[i] = p[2]
            uD[i] = u[2]
        elif(xI[i] < x4 and xI[i] > x3):
            rhoD[i] = rho[3]
            pD[i] = p[3]
            uD[i] = u[3]
        elif(xI[i] > x4):
            rhoD[i] = rho[4]
            pD[i] = p[4]
            uD[i] = u[4]
        eTemp = eqnState(pD[i],rhoD[i],uD[i],gamma)
        combList.append([pD[i],rhoD[i],uD[i],eTemp])
    analytSolFile(combList,fileName,fileStoreInfo)

def eqnState(p,rho,u,gamma):
    """Use this method to solve for pressure."""
    e = p/(gamma-1)+rho*u**2/2
    return e
def states(Gamma, gamma, beta, lBC, rBC,cL,cR):
    "Use this method to solve for the states"
    # Left side conditions
    pL = lBC[0]
    rhoL = lBC[1]
    uL = lBC[2]
    # Right side conditions
    pR = rBC[0]
    rhoR = rBC[1]
    uR = rBC[2]
    #iteratively solving
    tempList = pState3(Gamma,beta,lBC,rBC)
    #Pressure
    p2 = pL*(1-(gamma-1)/2*tempList[1]/cL)**(2*gamma/(gamma-1))
    p3 = tempList[0]
    p4 = p3
    #Velocity
    u2 = tempList[1]
    u4 = tempList[2]
    u3 = uR+(p3-pR)/np.sqrt(rhoR/2*((gamma+1)*p3+(gamma-1)*pR))

    #Density
    rho2 = rhoL*(1-(gamma-1)/2*u2/cL)**(2/(gamma-1))
    rho3 = rhoL*(p3/pL)**(1/gamma)
    rho4 = rhoR*(p4+Gamma*pR)/(pR+p4*Gamma)

    rho = [rhoL,rho2,rho3,rho4,rhoR]
    u = [uL,u2,u3,u4,uR]
    p = [pL,p2,p3,p4,pR]
    return [rho, u, p]

def pState3(Gamma, beta, lBC, rBC):
    """Use this method to iterate for pressure at state 3."""
    rP = [rBC[0],lBC[0]]
    bool =  True
    while(bool): #Loop to solve pressure iteratively
        pGuess = (rP[0]+rP[1])/2
        u4 = (pGuess-rBC[0])*np.sqrt((1-Gamma)/(rBC[1]*(pGuess+Gamma*rBC[0])))
        u2 = ((lBC[0]**beta-pGuess**beta)
            *np.sqrt((1-Gamma**2)*lBC[0]/(Gamma**2*lBC[1])))
        if(abs(u2-u4)<1e-15):
            bool = False
        elif(u2-u4 > 0):
            rP[0] = pGuess
        elif(u2-u4<0):
            rP[1] = pGuess

    return [pGuess,u2,u4]

def analytSolFile(data,fileName = None,extraInfo = None):
    """Use this method to write values to a file."""
    if(fileName is None):
        fileName = input("Enter file name and extension. ")
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

if __name__ == "__main__":
    eulerInfo = list()
    eulerInfo.append("Euler test case, time = 0")
    eulerInfo.append("dx = "+str(0.01))
    eulerInfo.append("dt = "+str(0.01))
    eulerInfo.append("nSteps = "+str(1))
    eulerInfo.append("[Pressure]   [Density]    [Velocity]    [Internal Energy]")
    t = np.linspace(0.1,.8,8)
    i = 0
    for time in t:
        str1 = "results/analyticalSol/test"+str(i)+".txt"
        i+=1
        sodShock(str1,eulerInfo,time)
