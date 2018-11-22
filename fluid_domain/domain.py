# Programmer: Anthony Walker
"""The purpose of this function is to generate a list of nodes"""
# READ ME:
# This file contains script for creating a multidimensional node tuple
# which acts as a domain object
# The class constuctor takes a fuile string or a multidimensional node tuple
# but not both. The domain ID is standardly set to zero.
# domain variables
# primaryDomain is a multidimensional node tuple containing the domains nodes
# lineDomain is a vector/flattened version of the primary domain.
# dims are the dimensions of the domain
# domainID is the id of the domain
# nRange is the range in a corresponding direction x, y, or z

from .node import node
import sys
import datetime as dt
import traceback
import numpy as np

class domain(object):  # A function to combine all methods
    """Use this class as an input node tuple for subdomain."""
#### Standard functions ###
    def __init__(self,fileStringIn=None,nodeTupleIn=None,domainID=0):
        """Use this as the class constructor."""
        #The purpose of this is to form a overloaded constructor of sorts
        try:
            if (isinstance(fileStringIn,str) and isinstance(nodeTupleIn,tuple)):
                print("In if")
                raise Warning("Warning:tuple ignored in domain generation")
            elif(isinstance(fileStringIn,str)):
                tempBool = True
            elif(isinstance(nodeTupleIn,tuple)):
                tempBool = False
            else:
                raise ValueError("ValueError: All values are of type None.")
        except Warning as tupIgWarn:
            print(tupIgWarn)
            tempBool = True
        except ValueError as vErr:
            print(vErr)
            sys.exit("Halting code execution.")
            tempBool = None
        finally:
            if tempBool is True:
                lineNodeTuple = self.NodeTupleGenerator(fileStringIn)
                self.primaryDomain = self.multidimNodeTuple(lineNodeTuple)
            elif tempBool is False:
                self.primaryDomain = nodeTupleIn
                lineNodeTuple = self.convMultiDim1D(self.primaryDomain)
            #Standard values for all constructors
            self.lineDomain = lineNodeTuple
            self.domainID = domainID
            self.dims = self.findTupleDimensions(lineNodeTuple)
            ranges = self.findTupleRange(lineNodeTuple)
            self.xRange = ranges[0]
            self.yRange = ranges[1]
            self.zRange = ranges[2]
    def __iadd__(self,other):
        """Use this method to add current values to domain values."""
        pass

    def __setitem__(self,key,item):
        """Use this method to copy domain or node values in assignment."""
        key = self.indexHandling(key)
        if isinstance(item,type(self)):
            self.copyAttributes(item)
        elif isinstance(item,type(np.array(1))):
            for i in range(key[0].start,key[0].stop+1):
                for j in range(key[1].start,key[1].stop+1):
                    self.primaryDomain[i][j].setValues(item)
        elif isinstance(item,type(self.primaryDomain[0][0])):
            for i in range(key[0].start,key[0].stop+1):
                for j in range(key[1].start,key[1].stop+1):
                    self.primaryDomain[i][j].setValues(item)

    def copyAttributes(self,domain):
        """Use this method to copy domain attributes."""
        self.lineDomain = domain.getLineDomain()
        self.domainID = domain.getDomainID()
        self.dims = domain.getDomainDims()
        ranges = domain.getDomainRanges()
        self.xRange = ranges[0]
        self.yRange = ranges[1]
        self.zRange = ranges[2]

    def __getitem__(self, index):
        """Use this method is to allow get indexing."""
        returnDomain = list()
        index = self.indexHandling(index)
        if index == False:
            return index
        for i in range(index[0].start,index[0].stop+1):
            for j in range(index[1].start,index[1].stop+1):
                #print(self.primaryDomain[i][j].getNodeLocation())
                returnDomain.append(self.primaryDomain[i][j])
        if len(returnDomain) != 1:
            returnDomain = self.multidimNodeTuple(tuple(returnDomain))
            returnDomain = domain(nodeTupleIn=returnDomain)
            return returnDomain
        else:
            return returnDomain[0]

    def indexHandling(self,index):
        """Use this method to handle indices for get and set item."""
        try:
            indexLength = len(index)

            if (indexLength == 2):
                dim = 2
            elif (indexLength == 3):
                dim = 3

            index = list(index)
            for i in range(dim):
                if not isinstance(index[i],slice):
                    index[i] = slice(index[i],index[i],1)
                #Statements to throw error for negative indices
                # if(index[i].start < 0 or index[i].stop < 0 ):
                #     raise IndexError("Negative index encountered")
                #     #sys.exit("Negative Index encountered, halting execution.")
        except TypeError as e:
            dim = 1
            if not isinstance(index,slice):
                index = slice(index,index,1)
            #Statements to throw error for negative indices
            # if(index.start < 0 or index.stop < 0 ):
            #     raise IndexError("Negative index encountered")
                #sys.exit("Negative Index encountered, halting execution.")
            index = [index]
        except IndexError as ie:
            return False

        for i in range(3-dim):
            index.append(slice(0,0,1))
        return index

#### Node handling functions ####

    def NodeTupleGenerator(self, fileString):
        """Use this function to generate a tuple from a txt file."""
        with open(fileString) as nodesFile:  # Opening file
            nodes = ()
            for i, C in enumerate(nodesFile):
                # Appending node tupcorners function,le
                nodeCoordStr = C.split(" ")
                nodeIDStr = nodeCoordStr[0]+"."+nodeCoordStr[1]+"."+nodeCoordStr[2]
                tempTup = (node(0, nodeIDStr, int(nodeCoordStr[0]), int(nodeCoordStr[1]),  int(nodeCoordStr[2])),)
                nodes = nodes+tempTup
            nodesFile.closed  # Closing file

        return nodes

    def multidimNodeTuple(self, nodeTuple):
        """Use this method to convert the node tuple generated to a 3d array."""

        dims = self.findTupleDimensions(nodeTuple)
        xTupList = list()  # Temporary list to convert to tuple
        #print(dims)
        for j in range(dims[0]):
            xTupList.append(nodeTuple[0:dims[1]])
            nodeTuple = nodeTuple[dims[1]:]
        xTup = tuple(xTupList)  # Creation of tuple from list.
        #print(xTup)
        return xTup

    def convMultiDim1D(self,tupToConv):
        """Use this method to convert a multidimensional node tuple to a 1D tuple."""
        returnTuple = list()
        for x in tupToConv:
            for y in x:
                returnTuple.append(y)
        return tuple(returnTuple)

    def findTupleDimensions(self, currNodeTuple):
        """Use this method to determine the dimensions of the tuple"""
        ranges = self.findTupleRange(currNodeTuple)

        xDim = abs(ranges[0].stop-ranges[0].start)+1
        yDim = abs(ranges[1].stop-ranges[1].start)+1
        zDim = abs(ranges[2].stop-ranges[2].start)+1
        self.xyEles = xDim*yDim
        self.xzEles = xDim*yDim
        self.yzEles = zDim*yDim
        return (xDim, yDim, zDim)

    def findTupleRange(self,currNodeTuple):
        """Use this method to determine the ranges of the tuple"""
        xTup = ()
        yTup = ()
        zTup = ()
        for k in currNodeTuple:  # Finding the max values to construct loops
            xTup = xTup + (k.getNodeLocation()[0],)
            yTup = yTup + (k.getNodeLocation()[1],)
            zTup = zTup + (k.getNodeLocation()[2],)
        xRange = range(min(xTup),max(xTup))
        yRange = range(min(yTup),max(yTup))
        zRange = range(min(zTup),max(zTup))
        return (xRange,yRange,zRange)

    def initializeDomain(self, *args):
        """Use this method to initialize the domain."""
        for arg in args:
            self.setNodeVals(arg[2],arg[0],arg[1])
        return

#### Verification Methods ####

    def printNodes(self,option = None):
        """Use this method to print all nodes for verification."""
        if option is not None:
            option = input("Would you like to print node values? [Y]es or [N]o. ")
        else:
            option = "Y"
        if option is "Y":
            for i in self.lineDomain:
                print(i.getNodeLocation(),i.getValues())
        return

    def printDomStats(self,option = None):
        if option is not None:
            option = input("Would you like to domain stats? [Y]es or [N]o. ")
        else:
            option = "Y"
        if option is "Y":
            print("*----------------*")
            print("Domain: "+str(self.domainID))
            print("Elements(xy): "+str(self.xyEles))
            #print("Elements(xz): "+str(self.xzEles))
            #print("Elements(yz): "+str(self.yzEles))
            print("----------------")
            print("Dimensions: "+str(self.dims))
            print("xRange: "+str(self.xRange))
            print("yRange: "+str(self.yRange))
            print("zRange: "+str(self.zRange))
            print("*----------------*")
        return

    def domainToFile(self,fileName = None,extraInfo = None):
        """Use this method to write node values to a file."""
        if(fileName is None):
            fileName = input("Enter file name and extensions. ")
        with open(fileName, 'w') as f:

            f.write("Date & Time: ")
            f.write(str(dt.datetime.now()))
            f.write("\n")

            if extraInfo is not None:
                for e in extraInfo:
                    f.write(e)
                    f.write("\n")
            for x in self:
                for y in x.getValues():
                    f.write(" %0.8f " % y)
                f.write("\n")
            f.closed

#### Get and Set Methods ####

    def setCorners():
        """Use this method to set the corners of the domain."""
        pass

    def setNodeVals(self, valueMat, xRange, yRange):
        """Use this function to set values of nodes"""
        # start is a tuple of the begining range Values
        # end is a tuple of the end range Values
        for i in xRange:
            for j in yRange:
                self.primaryDomain[i][j].setValues(valueMat)
        return

    def getPrimaryDomain(self):
        """Use this method to get the node domain."""
        return self.primaryDomain

    def setPrimaryDomain(self,domainToSet):
        """Use this method to set subdomain as a decomposed tuple."""
        self.primaryDomain = domainToSet
        return
    def getLineDomain(self):
        """Use this method to set subdomain as a decomposed tuple."""
        return self.lineDomain

    def getDomainDims(self):
        """Use this method to retrieve domain dimensions."""
        return self.dims

    def getDomainRanges(self):
        """Use this method to get domain ranges"""
        return (self.xRange,self.yRange,self.zRange)

    def setRange(self, x = None, y = None, z = None):
        """Use this method to set domain ranges"""
        if x is not None:
            self.xRange = x
        if y is not None:
            self.yRange = y
        if z is not None:
            self.zRange = z
        return

    def getDomainID(self):
        """Use this method to get domainID."""
        return self.domainID

    def setDomainID(self,idVal):
        """Use this method to set domainID."""
        self.domainID = idVal
        return

    def getElementSize(self):
        """Use this method to return the number of elements."""
        return (self.xyEles,self.xzEles,self.yzEles)
# Function main
if __name__ == "__main__":
    pass
