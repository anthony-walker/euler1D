# Programmer: Anthony Walker
"""The purpose of this file is to provide an applicable decomposition method."""

import numpy as np
# import time
from .domain import domain

def decompose(primaryDomain,numSubdomains):
    """Use this method to decompose a square."""
    dims = primaryDomain.getDomainDims()
    # m and n are the matrix dimensions
    m = dims[0]
    n = dims[1]

    finalPartDomain = tuple()

    if numSubdomains == 2:
        splitList = list()
        splitList.append(primaryDomain)
        splitDomains = domainSplit(splitList)
        finalPartDomain = finalPartDomain+(splitDomains[0][0],)+(splitDomains[0][1],)

    else:
        divs = domDivs(numSubdomains,m,n)

        domDiff = abs( numSubdomains-divs[0]*divs[1])
        if (m>n):
            partitionedDomain = divDecomp(m,divs[0],n,divs[1],primaryDomain,numSubdomains)
        else:
            partitionedDomain = divDecomp(m,divs[1],n,divs[0],primaryDomain,numSubdomains)
            i = 0


        partitionedDomain = tuple(partitionedDomain)
        
        if dimCheck(m,n):
            if domDiff>0:
                indices = findLargestDoms(domDiff,partitionedDomain)

                splitList = list();

                for i in indices:
                    splitList.append(partitionedDomain[i])

                splitDomains = domainSplit(splitList)

                indexCounter = 0
                for x in range(len(partitionedDomain)):
                    if x == indices[indexCounter]:
                        #print("Adding two")
                        finalPartDomain = finalPartDomain+splitDomains[indexCounter]
                        indexCounter = indexCounter+1
                    else:
                        finalPartDomain = finalPartDomain+(partitionedDomain[x],)
            else:
                finalPartDomain = partitionedDomain
        else:
            finalPartDomain = partitionedDomain
    k = 0
    for x in finalPartDomain:
        x.setDomainID(k)
        k = k+1
    return finalPartDomain

def divDecomp(m, div1, n, div2, primaryDomain, numSubdomains):
    """Use this method to determine 1D or 2D."""
    if not dimCheck(m,n):
        pDom = divDecomp1D(m, div1,primaryDomain, numSubdomains)
    else:
        pDom = divDecomp2D(m, div1, n, div2, primaryDomain, numSubdomains)
    return pDom

def dimCheck(m,n):
    if m == 1 or n == 1:
        bool = False
    else:
        bool = True
    return bool

def divDecomp2D(m, div1, n, div2, primaryDomain, numSubdomains):
    """Use this method to decompose the domain by div1 and div2"""
    mDiv = int(np.floor(m/div1))
    nDiv = int(np.floor(n/div2))
    ranges = primaryDomain.getDomainRanges()
    xMax = ranges[0].stop
    yMax = ranges[1].stop
    partitionedDomain = list()

    for x in range(div1-1):
        for y in range(div2-1):
                #print(x,y)
                tempDom = primaryDomain[mDiv*(x):mDiv*(x+1)-1,nDiv*(y):nDiv*(y+1)-1]
                partitionedDomain.append(tempDom)
                #tempDom.printDomStats()

                if (y == div2-2):
                    #print(y+1,div2,"loop1")
                    tempDom = primaryDomain[mDiv*(x):mDiv*(x+1)-1,nDiv*(y+1):yMax]
                    partitionedDomain.append(tempDom)
                    #tempDom.printDomStats()
    for y in range(div2-1):
            tempDom = primaryDomain[mDiv*(div1-1):xMax,nDiv*(y):nDiv*(y+1)-1]
            #print(x+1,div1,"loop2")
            partitionedDomain.append(tempDom)
            #tempDom.printDomStats()
            if (y+1 == div2-1):
                tempDom = primaryDomain[mDiv*(div1-1):xMax,nDiv*(div2-1):yMax]
                partitionedDomain.append(tempDom)
                #tempDom.printDomStats()
    return partitionedDomain

def divDecomp1D(m,div1, primaryDomain, numSubdomains):
    mDiv = int(np.floor(m/numSubdomains))
    ranges = primaryDomain.getDomainRanges()
    xMax = ranges[0].stop
    partitionedDomain = list()

    for x in range(numSubdomains-1):
        #print(x,y)
        tempDom = primaryDomain[mDiv*(x):mDiv*(x+1)-1,0]
        partitionedDomain.append(tempDom)
        #tempDom.printDomStats()
        if(x+1 == numSubdomains-1):
            tempDom = primaryDomain[mDiv*(numSubdomains-1):xMax,0]
            partitionedDomain.append(tempDom)
            #tempDom.printDomStats()

    return partitionedDomain
def domDivs(numSubdomains,m,n):
    """Use this method to get the domain divisions."""

    rootOfDoms = np.sqrt(numSubdomains)
    if rootOfDoms > np.floor(rootOfDoms)+1/2:

        div1 = int(np.floor(rootOfDoms)+1)
        div2 = int(np.floor(rootOfDoms))
        # div1 = divMin(div1,m)
        # div2 = divMin(div2,n)
    else:
        div1 = int(np.floor(rootOfDoms))
        div2 = int(np.floor(rootOfDoms))
        # div1 = divMin(div1,m)
        # div2 = divMin(div2,n)
    return (div1,div2)

def divMin(currDiv,dim):
    "Use this method to return the appropriate div value"
    if(currDiv>dim):
        retVal = dim
    else:
        retVal = currDiv
    return retVal

def findLargestDoms(domDiff,pDomain):
    """Use this method to find the largest domains."""

    indexList = list(range(0,len(pDomain)))
    elementList = list()
    returnIndices = list()
    for x in pDomain:
        tempEleSize = x.getElementSize()
        tempEleSize = tempEleSize[0]
        elementList.append(tempEleSize)

    domainList = zip(elementList,indexList)
    domainList = sorted(domainList)
    domainList = reversed(domainList)
    i = 0
    for x in domainList:
        if i < domDiff:
            tempVal = x[1]
            returnIndices.append(tempVal)
        i = i+1
    returnIndices.sort();
    return returnIndices

def domainSplit(splitList,splitDir = 1):
    """Use this Doms:thod to split a domain."""
    returnList = list()

    if splitDir == 1:
        for dom in splitList:
            dims = dom.getDomainDims()
            dom1 = dom[0:dims[0]-1,0:int(np.floor(dims[1]/2))-1]
            dom2 = dom[0:dims[0]-1,int(np.floor(dims[1]/2)):dims[1]-1]
            returnList.append((dom1,dom2))
    elif splitDir == 2:
        for dom in splitList:
            dims = dom.getDomainDims()
            dom1 = dom[0:int(np.floor(dims[0]/2))-1,0:dims[1]-1]
            dom2 = dom[int(np.floor(dims[0]/2)):dims[0]-1,0:dims[1]-1]
            returnList.append((dom1,dom2))

    return tuple(returnList)

if __name__ == "__main__":
    testDomain = domain("nodes.txt")
    testDomain = decompose(4,testDomain)
