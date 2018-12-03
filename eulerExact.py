#Programmer: Anthony Walker
import numpy as np
import datetime as dt
#Allocation of global variables
gamma = 0
alpha = 0
beta = 0
cL = 0
cR = 0
lbc = 0
rbc = 0
t = 0

def eulerExact(fileName,fileStoreInfo,time = 0.2,lBC = (1.0,1.0,0,2.5),rBC = (0.1,0.125,0,0.25),g = 1.4,numPts = 1000):
    """Use this method to determine the analytic solution of sod shock problem."""
    #Important coeffcients
    global gamma
    gamma = g
    global alpha
    alpha = (gamma+1)/(gamma-1)
    global beta
    beta = (2*gamma)/(gamma-1)
    #Boundary conditions
    global lbc
    lbc = lBC
    global rbc
    rbc = rBC
    #Time
    global t
    t = time
    #Speed of sound for states 1 and 5
    global cL
    cL = np.sqrt(gamma*lbc[0]/lbc[1])
    global cR
    cR = np.sqrt(gamma*rbc[0]/rbc[1])
    #Determining states
    statesList = states2and3()
    rho = statesList[0]
    u = statesList[1]
    p = statesList[2]
    speed = statesList[3]
    #Location
    x1 = 0
    xshock = x1+speed[3]*t
    xcontact = x1 +speed[2]*t
    xfR = x1 + speed[1]*t
    xfL = x1 + speed[0]*t
    x = [xfL,xfR,x1,xcontact,xshock]

    #Getting Data
    xI = np.linspace(-0.5,0.5,numPts)
    pD = np.zeros(numPts)
    rhoD = np.zeros(numPts)
    uD = np.zeros(numPts)
    combList = list()
    for i in range(numPts): #State L
        if(xI[i]<xfL):
            pD[i] = lbc[0]
            rhoD[i] = lbc[1]
            uD[i] = lbc[2]
        elif(xI[i] < xfR and xI[i] > xfL): #State 4
            pD[i] = lbc[0]*(1+(xfL-xI[i])/(cL*alpha*t))**(beta)
            rhoD[i] = lbc[1]*(1+(xfL-xI[i])/(cL*alpha*t))**(beta)
            uD[i] = lbc[2]-2/(gamma+1)*(xfL-xI[i])/t
        elif(xI[i] < xcontact and xI[i] > xfR): #State 3
            rhoD[i] = rho[1]
            pD[i] = p[1]
            uD[i] = u[1]
        elif(xI[i] < xshock and xI[i] > xcontact): #State 2
            rhoD[i] = rho[0]
            pD[i] = p[0]
            uD[i] = u[0]
        elif(xI[i] > xshock):
            pD[i] = rbc[0]
            rhoD[i] = rbc[1]
            uD[i] = rbc[2]
        eTemp = eqnState(pD[i],rhoD[i],uD[i],gamma)
        combList.append([pD[i],rhoD[i],uD[i],eTemp])
    analytSolFile(combList,fileName,fileStoreInfo)

def eqnState(p,rho,u,gamma):
    """Use this method to solve for pressure."""
    e = p/(gamma-1)+rho*u**2/2
    return e

def states2and3(tol = (10**-10)):
    "Use this method to solve for the states"
    # Left side conditions
    pL = lbc[0]
    rhoL = lbc[1]
    uL = lbc[2]
    # Right side conditions
    pR = rbc[0]
    rhoR = rbc[1]
    uR = rbc[2]
    #iteratively solving Rankine-Hugoniot relation
    P = RankineHugoniot(2,pL,pR,tol)
    #State 2
    rho2 = (1+beta*P)/(beta+P)*rhoR
    p2 = P*pR
    u2 = uL+2*cL/(gamma-1)*(1-(P*pR/pL)**(1/beta))
    ss = uR+cR*(((gamma-1)+(gamma+1)*P)/(2*gamma))*(1/2) #Shock loc

    #State 3
    p3 = p2
    u3 = u2
    rho3 = rhoL*(P*(pR/pL))**(1/gamma)
    sc = uL+2*cL/(gamma-1)*(1-(P*pR/pL)**(1/beta)) #Contact loc

    #State 4
    sfL = -cL #Expansion fan loc left
    sfR = u2-cL #Expansion fan loc right

    #Return variables
    loc =[sfL,sfR,sc,ss]
    rho = [rho2,rho3]
    u = [u2,u3]
    p = [p2,p3]
    return [rho, u, p,loc]

def RankineHugoniot(P,pL,pR,tol):
    """Use this method to solve the Rankine-Hugoniot relation."""
    prevP = P+1
    while(abs(P-prevP)>tol):
        prevP = P
        P = (pL/pR*(1-(gamma-1)*(cR/cL)*(P+1)/
        np.sqrt(2*gamma*(2*gamma+(gamma+1)*(P-1))))**beta)
    return P

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
    str1 = "results/aSol/test1.txt"
    # t = np.linspace(0.1,.8,8)
    # i = 0
    # for time in t:
    #     str1 = "results/analyticalSol/test"+str(i)+".txt"
    #     i+=1
    eulerExact(str1,eulerInfo,time = 0.2)
