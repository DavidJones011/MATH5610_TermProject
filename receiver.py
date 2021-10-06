from MATH5610_TermProject.satellite import magnitude, satelliteTimeAndLocOnSend
import sys
import math
import os
import re
from decimal import *

def main() :

    global pi, s, c, r
    pi, s, c, r = getData()
    sattellites, startIndicies = readSatelliteData()
    num = len(sattellites)

    for i in range(0,num) : 

        x_v = [0,0,0]
        count = 0
        if (i + 1) < num :
            count = (startIndicies[i+1] - startIndicies[i]) + 1
        else :
            count = (num - startIndicies[i]) + 1

        calculateFirstOrderPartDeriv(x_v, sattellites, startIndicies[i], count)
        #calculateSecOrderPartDeriv(x_v, satellites, startIndicies[i], count)

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
            pi = Decimal(tokens[0])
        elif(tokens[2] == 's') :
            s = Decimal(tokens[0])
        elif(tokens[2] == 'c') :
            c = Decimal(tokens[0])
        elif(tokens[2] == 'R') :
            r = Decimal(tokens[0])

    data_file.close()
    return pi, s, c, r

# reads the satellite data from std::in
def readSatelliteData() :
    satellites = list()
    startIndices = list()

    prevIndex = 0
    curIndex = 0

    index = 0
    for line in sys.stdin :
        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        cur_satellite = [0, 0, [0,0,0]]
        cur_satellite[0] = int(tokens[0])
        cur_satellite[1] = Decimal(tokens[1])
        cur_satellite[2][0] = Decimal(tokens[2])
        cur_satellite[2][1] = Decimal(tokens[3])
        cur_satellite[2][2] = Decimal(tokens[4])
        satellites.append(cur_satellite)

        if cur_satellite[0] < prevIndex : 
            startIndices.append(index)

        prevIndex = cur_satellite[0]
        index = index + 1

    return satellites, startIndices

def calculateFirstOrderPartDeriv(vec, satellites, start, count) :

    for i in range(start, end+1) :


    c_magnitude = magnitude(subVectors(c_satellite[2], vec))
    n_magnitude = magnitude(subVectors(n_satellite[2], vec))

    a_i = n_magnitude - c_magnitude - c * (c_satellite[1] - n_satellite[1])

    der = addVectors(invScaleVector(-n_magnitude, subVectors(n_satellite[2], vec)),
                    invScaleVector(c_magnitude, subVectors(c_satellite[2], vec)))

    der = scaleVector(a_i, der)

    return der

def calculateSecOrderPartDeriv(vec, satellites, start, count) :
    pass

# returns the magnitude of u
def magnitude(u):
    return Decimal.sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]))

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
    x = Decimal(vec[0] * Decimal(scale))
    y = Decimal(vec[1] * Decimal(scale))
    z = Decimal(vec[2] * Decimal(scale))
    return [x, y, z]

# returns a vector that is scaled by the inverse of the scalar
def invScaleVector(scale, vec) :
    x = Decimal(vec[0] / Decimal(scale))
    y = Decimal(vec[1] / Decimal(scale))
    z = Decimal(vec[2] / Decimal(scale))
    return [x, y, z]

# creating the Jacobian J(x) as defined in exercise 13
# we will use J(x) to solve the nonlinear system of 4 equations
# there is probably an elegant way to automate the creation of the jacobian and append each element via a for loop

J = [[(xS1 - x)/magnitude(xS1 - x) - (xS2 - x)/magnitude(xS2 - x), (yS1 - y)/magnitude(xS1 - x) - (yS2 - y)/magnitude(xS2 - x), (zS1 - z)/magnitude(xS1 - x) - (zS2 - z)/magnitude(xS2 - x)],
     [(xS2 - x)/magnitude(xS2 - x) - (xS3 - x)/magnitude(xS3 - x), (yS2 - y)/magitude(xS2 - x) - (yS3 - y)/magnitude(xS3 - x), (zS2 - z)/magnitude(xS2 - x) - (zS3 - z)/magnitude(xS3 - x)],
     [(xS3 - x)/magnitude(xS3 - x) - (xS4 -x)/magnitude(xS4 - x), (yS3 - y)/magnitude(xS3 - x) - (yS4 - y)/magnitude(xS4 - x), (z3 - z)/magnitude(xS3 - x) - (z4 - z)/magnitude(xS4 - x)]]

# we want to solve J(x^(k))s^(k) = -F(x^(k)) where x^(k+1) = x^(k) + s^(k)
#
#
#
#

# creating derivatives as outlined in exercise 14
# going to us xSi, can we refactor as needed and fill in for i such that i = 1,2,3,4...m-1?
# iPlus1 --> i + 1, when referring to indices, not sure how we want assign values here, a function with a loop may be better?

# initialize weird variables, then can loop through and update them accordinly

# using a null initial value
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

# def DYiDz(NiPlus1, ySiPlus1, zSiPlus1, Ni, ySi, zSi, y, z):

# def DZiDz(NiPlus1, zSiPlus1, Ni, zSi, z):
    
main()
