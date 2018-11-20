# Programmer: Anthony Walker
"""The purpose of this class is to act as a node point in the subdomain."""
# READ ME:
# This file contains the script for the class "Node"
# The node class is used to represent discretized points in a fluid domain
# which is represented by a domain object.
# node variables
# Values is a list of values corresponding to a Node
# pressure, density, temperature, velocity, etc..
# nodeID is the nodes identifier

import numpy as np

class node(object):
    """Use the node class to represent points in the subdomain."""

#### Standard Functions ####
    def __init__(self, values, nodeID, x, y, z):  # Node class constructor
        """Use this definition as a class constructor."""
        self.values = np.array(values)
        self.nodeID = nodeID
        self.loc = (x, y, z)
    def __iadd__(self,other):
        """Use this method to add current values to domain values."""
        self.values += other
        return self
    def __setitem__(self,key,values):
        self.values[key] = values
    def __getitem__(self,key):
        return self.values[key]

#### Get and Set functions ####
    def getNodeID(self):  # accessor for nodeID
        """Use this method in order to get nodeID."""
        return self.nodeID

    def getValues(self):  # accessor for values
        """Use this method to get values stored in the node."""
        return self.values

    def setNodeID(self, value):  # mutator for nodeID
        """Use this method to set nodeID of the node."""
        self.nodeID = value
        return

    def setValues(self, values):
        """Use this method to set values stored in the node."""
        self.values = np.array(values)
        return

    def getNodeLocation(self):
        """Use this method to get the node location."""
        return self.loc

    def setNodeLocation(self, x, y, z):
        """Use this method to set node location."""
        self.loc = (x, y, z)
        return

#### other functions ####
    def isBoundary(self):
        """Use this method to check if the node is on a boundary."""
        return self.boundBool

if __name__ == "__main__":
    pass
