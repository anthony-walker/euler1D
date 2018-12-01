# Programmer: Anthony Walker
"""The purpose of this file is post processing of data obtain
from the fluid domain."""
#READ ME:
#Use this to plot all of the data in a folder that is created by the
#euler1D code.

import matplotlib as mpl
mpl.use("Tkagg")
import matplotlib.pyplot as plt
import re
import numpy as np
import os
def dataPlot(directory,labelList,l2s):
    """Use this method to plot data appropriately."""
    fileList = sorted(os.listdir(directory))
    fileDataList = list()
    bool = True
    for i in fileList:
        columnList = convFileToPlotData(directory+i,l2s)
        if bool:
            for x in range(len(columnList)):
                fileDataList.append(tuple())
            bool = False
        for j in range(len(columnList)):
            fileDataList[j]+=(columnList[j],)
    x =np.linspace(0,1,len(columnList[0]))
    fig, axes = plt.subplots(2,2)
    i = 0
    axes[0,0].set_title(labelList[0])
    axes[0,1].set_title(labelList[1])
    axes[1,0].set_title(labelList[2])
    axes[1,1].set_title(labelList[3])

    for y in fileDataList[0]:
        axes[0,0].plot(x,y)
    for y in fileDataList[1]:
        axes[0,1].plot(x,y)
    for y in fileDataList[2]:
        axes[1,0].plot(x,y)
    for y in fileDataList[3]:
        axes[1,1].plot(x,y)
    #axes[0,0].plt.plot(x,y)
    # figure manager for Tkagg

    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    # plt.ylabel("variable")
    # plt.xlabel("x")
    plt.show()

def convFileToPlotData(fileName,linesToSkip):
    """"Use this method to convert domain file to plot data."""

    columnList = list()
    with open(fileName, 'r') as f:
        for x in range(linesToSkip):
            f.readline()
        currLine = f.readline()
        currLine = [float(s) for s in re.findall(r'-?\d+\.?\d*', currLine)]
        for n in range(len(currLine)):
            columnList.append(list())
            columnList[n].append(currLine[n])
        for i in f:
            i = [float(s) for s in re.findall(r'-?\d+\.?\d*', i)]
            for n in range(len(i)):
                columnList[n].append(i[n])
        f.closed
        return columnList
if __name__ == "__main__":
    labels = ["Pressure","Density","Velocity","Internal Energy"]
    lSkip = 6
    dataPlot("results/testing/",labels,lSkip)
