from MATH5610_TermProject.satellite import magnitude, satelliteTimeAndLocOnSend
import sys
import mpmath as mp
import os
import re

def main() :
    # set the decimal place for calculations
    mp.mp.dps = 17
    # create and set global variables
    global pi, s, c, r
    pi, s, c, r = getData()
    sattellites, startIndicies = readSatelliteData()
    numPositions = len(startIndicies) - 1

    for i in range(0, numPositions) :

        startIndex = startIndicies[i]
        count = startIndicies[i+1] - startIndex

        for i in range(0, count) :

            # calculate first order and second order partial derivatives
            

            # newtons method

            pass

        # convert to geodesic coordinates

        # ouput results

        pass

    pass

# grab data from data.dat in local folder
def getData() :
    # initialize satellite data
    satellites = list()
    for i in range(0, 24):
        satellites.append([[0,0,0],[0,0,0],0,0,0])

    data_file = open(os.path.join(sys.path[0], 'data.dat'), "r")

    for i in range(0,4) :

        line = data_file.readline()

        if line == '':
            break

        tokens = re.split('\s|,', line)

        # remove empty strings from the list of tokens
        while("" in tokens) :
            tokens.remove("")

        if(len(tokens) == 0) :
            break

        if(tokens[2] == 'pi') :
            pi = mp.mpf(tokens[0])
        elif(tokens[2] == 's') :
            s = mp.mpf(tokens[0])
        elif(tokens[2] == 'c') :
            c = mp.mpf(tokens[0])
        elif(tokens[2] == 'R') :
            r = mp.mpf(tokens[0])

    data_file.close()
    return pi, s, c, r

# reads the satellite data from std::in
def readSatelliteData() :
    satellites = list()
    startIndices = list()
    #satellite indices range from 0-23
    prevIndex = 24
    index = 0

    for line in sys.stdin :
        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        cur_satellite = [0, 0, [0,0,0]]
        cur_satellite[0] = int(tokens[0])
        cur_satellite[1] = mp.mpf(tokens[1])
        cur_satellite[2][0] = mp.mpf(tokens[2])
        cur_satellite[2][1] = mp.mpf(tokens[3])
        cur_satellite[2][2] = mp.mpf(tokens[4])
        satellites.append(cur_satellite)

        # found a new set of satellites for a new position
        # store the starting index
        if cur_satellite[0] < prevIndex : 
            startIndices.append(index)

        prevIndex = cur_satellite[0]
        index = index + 1

    startIndices.append(len(startIndices))
    return satellites, startIndices

def calculateFirstOrderPartDeriv(vec, satellites, start, count) :

    

    pass

def calculateSecOrderPartDeriv(vec, satellites, start, count) :
    pass

# returns the magnitude of u
def magnitude(u):
    return mp.sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]))

# returns the dot product of two vectors u and v
def dotProduct(u, v) :
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2])

# adds two vectors
def addVectors(u,v) :
    return [u[0] + v[0], u[1] + v[1], u[2] + v[2]]

# subtracts two vectors
def subVectors(u,v) :
    return [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

# returns a vector that is scaled by the scalar
def scaleVector(scale, vec) :
    x = mp.mpf(vec[0] * mp.mpf(scale))
    y = mp.mpf(vec[1] * mp.mpf(scale))
    z = mp.mpf(vec[2] * mp.mpf(scale))
    return [x, y, z]

# returns a vector that is scaled by the inverse of the scalar
def invScaleVector(scale, vec) :
    x = mp.mpf(vec[0] / mp.mpf(scale))
    y = mp.mpf(vec[1] / mp.mpf(scale))
    z = mp.mpf(vec[2] / mp.mpf(scale))
    return [x, y, z]

# creating the Jacobian J(x) as defined in exercise 13
# we will use J(x) to solve the nonlinear system of 4 equations
# there is probably an elegant way to automate the creation of the jacobian and append each element via a for loop

#J = [[(xS1 - x)/magnitude(xS1 - x) - (xS2 - x)/magnitude(xS2 - x), (yS1 - y)/magnitude(xS1 - x) - (yS2 - y)/magnitude(xS2 - x), (zS1 - z)/magnitude(xS1 - x) - (zS2 - z)/magnitude(xS2 - x)],
#     [(xS2 - x)/magnitude(xS2 - x) - (xS3 - x)/magnitude(xS3 - x), (yS2 - y)/magitude(xS2 - x) - (yS3 - y)/magnitude(xS3 - x), (zS2 - z)/magnitude(xS2 - x) - (zS3 - z)/magnitude(xS3 - x)],
#     [(xS3 - x)/magnitude(xS3 - x) - (xS4 -x)/magnitude(xS4 - x), (yS3 - y)/magnitude(xS3 - x) - (yS4 - y)/magnitude(xS4 - x), (z3 - z)/magnitude(xS3 - x) - (z4 - z)/magnitude(xS4 - x)]]

# we want to solve J(x^(k))s^(k) = -F(x^(k)) where x^(k+1) = x^(k) + s^(k)
#
#
#
#

# creating derivatives as outlined in exercise 14
# going to us xSi, can we refactor as needed and fill in for i such that i = 1,2,3,4...m-1?
# iPlus1 --> i + 1, when referring to indices, not sure how we want assign values here, a function with a loop may be better?

# initialize weird variables, then can loop through and update them accordinly
xSi = None
xSiPlus1 = None
ySi = None
ySiPlus1 = None
zSi = None
zSiPlus1 = None
tSi = None
tSiPlus1 = None
NiPlus1 = None

# we make some generalizations here to make the first order partial derivatives easier to write
Ni = magnitude(xSi - x)
Ai = NiPlus1 - Ni - c(tSi - tSiPlus1)
Xi = -(xSiPlus1 - x)/(NiPlus1) + (xSi - x)/Ni
Yi = -(ySiPlus1 - y)/(NiPlus1) + (ySi - y)/Ni
Zi = -(zSiPlus1 - z)/(NiPlus1) + (zSi - z)/Ni

# just writing the first order partials as seen in equation (73) on hw01
def DXiDx(NiPlus1, xSiPlus1, Ni, xSi, x):
    return ((NiPlus1^2) - (xSiPlus1 - x)^2)/(NiPlus1^3) - ((Ni^2) - (xSi - x)^2)/(Ni^3) 

def DYiDx(NiPlus1, ySiPlus1, xSiPlus1, Ni, ySi, xSi, x, y):
    return -((ySiPlus1 - y) * (xSiPlus1 - x))/(NiPlus1^3) + ((ySi - y) * (xSi - x))/(Ni^3) #equivalent to DXiDy (partial of Xi with respect to y)

def DXiDz(NiPlus1, xSiPlus1, zSiPlus1, Ni, xSi, zSi, x, z):
    return -((xSiPlus1 - x) * (zSiPlus1 - z))/(NiPlus1^3) + ((xSi - x) * (zSi - z))/(Ni^3) #equivalent to DZiDx

def DYiDy(NiPlus1, ySiPlus1, Ni, ySi, y):
    return ((NiPlus1^2) - (ySiPlus1 - y)^2)/(NiPlus1^3) - ((Ni^2) - (ySi - y)^2)/(Ni^3)

def DYiDz(NiPlus1, ySiPlus1, zSiPlus1, Ni, ySi, zSi, y, z):
    return -((ySiPlus1 - y) * (zSiPlus1 - z))/(NiPlus1^3) + ((ySi - y) * (zSi - z))/(Ni^3)) #equivalent to DZiDy

def DZiDz(NiPlus1, zSiPlus1, Ni, zSi, z):
    return ((NiPlus1^2) - (zSiPLus1 - z)^2)/(NiPlus1^3) - ((Ni^2) - (zSi - z)^2)/(Ni^3)
    
main()
