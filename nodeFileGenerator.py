# Programmer: Anthony Walker
"""The purpose of this is to generate rectangular node files."""
#READ ME:
#This can generate text files to be used in domain creation.
#It was used to generate the example euler1D for domain creation.

def generateNodeFile(fileName, xRange=range(0,1), yRange=range(0,1), zRange = range(0,1)):

    with open(fileName,"w") as nodesFile:

        for i in xRange:
            for j in yRange:
                for k in zRange:
                    writeStr = str(i)+" "+str(j)+" "+str(k)+"\n"
                nodesFile.write(writeStr)
        nodesFile.closed


if __name__ == "__main__":

    generateNodeFile("euler1D.txt", range(0,101), range(0,1))
