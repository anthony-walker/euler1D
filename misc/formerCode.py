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
