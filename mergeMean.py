#-------------------------------------------------------------------------------
# Name:        mergeMean.py
# Purpose:     Smooth probe beta values within windows
#
# Author:      Dan
#
# Created:     27/12/2013
# Copyright:   (c) Dan 2013
#-------------------------------------------------------------------------------

def main():
    myFile = "no_quote.txt"
    numRows = 0
    numCols = 0
    with open(myFile,"r") as cntStuff:
        for line in cntStuff:
            if(numRows == 0):
                header = line.split()
                numCols = len(header)
            numRows = numRows + 1
    matrix = [[0.0 for x in xrange(numCols-1)] for x in xrange(numRows-1)]
    genes = ["NA" for x in xrange(numRows-1)]
    curRow = 0
    with open(myFile,"r") as loadVals:
        for line in loadVals:
            if(curRow > 0) :
                tempSplit = line.split()
                for x in xrange(numCols):
                    if x < 3:
                        matrix[curRow-1][x] = int(tempSplit[x])
                    elif x > 3:
                        matrix[curRow-1][x-1] = float(tempSplit[x])
                    else:
                        genes[curRow-1] = tempSplit[x]
            curRow = curRow + 1
    dupMat = [[0.0 for x in xrange(numCols-4)] for x in xrange(numRows-1)]
    baseRow = 0
    capRow = 0
    trim = .1
    dist = 500
    for row in xrange(len(dupMat)):
        while ((matrix[baseRow][2] + dist < matrix[row][2]) or (matrix[baseRow][1] != matrix[row][1])):
            baseRow = baseRow + 1
	capRow = row
        while (matrix[capRow][2] < matrix[row][2] + dist) and (capRow+1 < len(dupMat)):
            capRow = capRow + 1
        while (matrix[capRow][1] != matrix[row][1]) or (matrix[capRow][2] > matrix[row][2]+dist):
            capRow = capRow - 1
        indexTrim = int((capRow-baseRow + 1) * trim)
        baseRow = baseRow + indexTrim
        capRow = capRow - indexTrim
        divisor = 0.0
        for i in xrange(baseRow,capRow+1):
            divisor = divisor + 1.0
            for j in xrange(len(dupMat[0])):
                dupMat[row][j] = dupMat[row][j] + matrix[i][j+3]
        for j in xrange(len(dupMat[0])):
	    if(divisor == 0):
		print(baseRow,capRow,indexTrim,row)
            dupMat[row][j] = dupMat[row][j] / divisor
        baseRow = baseRow - indexTrim
        capRow = capRow + indexTrim
    outputfile = "RReadyFinal.txt"
    outStream = open(outputfile,"w")
    nextLine = "\"" + header[0] + "\""
    for i in xrange(1,len(header)):
        nextLine = nextLine + " \"" + header[i] + "\""
    nextLine = nextLine + "\n"
    outStream.write(nextLine)
    for row in xrange(len(matrix)):
        nextLine = str(matrix[row][0]) + " " + str(matrix[row][1]) + " " + str(matrix[row][2]) + " \"" + genes[row] + "\""
        for col in xrange(len(dupMat[row])):
            nextLine = nextLine + " " + ("%.4f" % dupMat[row][col])
        outStream.write(nextLine +"\n")
    outStream.close()









if __name__ == '__main__':
    main()
