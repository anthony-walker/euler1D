#Old minmod method from euler1D
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
#From euler 1D
if(eE):
    eulEStr = dirStr+"/eESol1"+eStr+".txt"
    euE.eulerExact(eulEStr,eulerInfo,tCurr,numPts = dims[0])
if(sA):
    sodStr = dirStr+"/aSol1"+eStr+".txt"
    sd.sodShock(sodStr,eulerInfo,tCurr,numPts = dims[0])

#From euler1D
def RK2(fcn,domain):
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
    uL = qL[1]/qL[0]
    uR = qR[1]/qR[0]
    eL = qL[2]/qL[0]
    eR = qR[2]/qR[0]
    rootrhoL = np.sqrt(qL[0])
    rootrhoR = np.sqrt(qR[0])
    qSP[0] = rootrhoL*rootrhoR
    denom = 1/(rootrhoL+rootrhoR)
    qSP[1] = (rootrhoL*uL+rootrhoR*uR)*denom
    qSP[2] = (rootrhoL*eL+rootrhoR*eR)*denom
    pSP = eqnState(qSP[0],qSP[1],qSP[2])
    rSP = np.sqrt(gamma*pSP/qSP[0])+abs(qSP[1])
    return rSP*(qL-qR)

def eqnState(rho,u,e):
    """Use this method to solve for pressure."""
    P = (gamma-1)*(rho*e-rho*u*u/2)
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
    if len(nV) == 1:
        nV = nV[0]

    return nV

def makeQ(nodes):
    """Use this method to obtain Q values from nodeValues."""
    Q = tuple()
    #P, rho, rho*u, rho*e
    for x in nodes:
        Q += (np.array([x[0],x[1],x[1]*x[2],x[1]*x[3]]),)
    if len(Q) == 1:
        Q = Q[0]
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
