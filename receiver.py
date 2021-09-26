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

main()