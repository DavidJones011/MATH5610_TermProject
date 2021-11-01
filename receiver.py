import sys
import mpmath as mp
import numpy as np
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
    max_iterations = 200

    output_file = open(os.path.join(sys.path[0], 'satellite.log'), "w")

    for i in range(0, numPositions) :

        startIndex = startIndicies[i]
        count = startIndicies[i+1] - startIndex
        x_i = np.array([-1795225.28, -4477174.36, 4158593.45])

        for k in range (0, max_iterations) :   
            F, J_inv = test(x_i, sattellites, startIndex, count)
            s_i = np.linalg.solve(J_inv, -F)
            
            x_ip1 = x_i + s_i
            diff = x_ip1 - x_i
            x_i = x_ip1
            if(abs(diff) < 0.01) :
                break
            pass

        sys.stdout.write("{} {} {}\n".format(x_i[0], x_i[1], x_i[2]))
        # convert to geodesic coords and print
        #output_file.write("{} {} {} {} {}\n".format(i, t_s, x_s[0], x_s[1], x_s[2]))
        #sys.stdout.write("{} {} {} {} {}\n".format(i, t_s, x_s[0], x_s[1], x_s[2]))

        pass

    output_file.close()
    pass

# grab data from data.dat in local folder
def getData() :
    # initialize satellite data
    #satellites = list()
    #for i in range(0, 24):
    #    satellites.append([[0,0,0],[0,0,0],0,0,0])

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
            pi = np.float64(tokens[0])
        elif(tokens[2] == 's') :
            s = np.float64(tokens[0])
        elif(tokens[2] == 'c') :
            c = np.float64(tokens[0])
        elif(tokens[2] == 'R') :
            r = np.float64(tokens[0])

    data_file.close()
    return pi, s, c, r

# reads the satellite data from std::in
def readSatelliteData() :
    satellites = list()
    startIndices = list()
    prevIndex = 24
    index = 0

    for line in sys.stdin :
        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        cur_satellite = [0, 0, np.array([0,0,0])]
        cur_satellite[0] = int(tokens[0])
        cur_satellite[1] = np.float64(tokens[1])
        cur_satellite[2] = np.array([np.float64(tokens[2]), np.float64(tokens[3]), np.float64(tokens[4])])
        #cur_satellite[2][0] = mp.mpf(tokens[2])
        #cur_satellite[2][1] = mp.mpf(tokens[3])
        #cur_satellite[2][2] = mp.mpf(tokens[4])
        satellites.append(cur_satellite)

        # found a new set of satellites for a new position
        # store the starting index
        if cur_satellite[0] <= prevIndex : 
            startIndices.append(index)

        prevIndex = cur_satellite[0]
        index = index + 1

    startIndices.append(len(satellites))
    return satellites, startIndices

#
def test(x, satellites, startIndex, count) :
    
    #initialize local variables
    DFdx = 0.0
    DFdy = 0.0
    DFdz = 0.0

    J_inv = np.array([[0, 0, 0],
                      [0, 0, 0],
                      [0, 0, 0]])

    for i in range(0, count-1) :
        diff = satellites[i][2] - x
        diff2 = satellites[i+1][2] - x

        N_i = np.linalg.norm(diff)
        N_ip1 = np.linalg.norm(diff2)
        Ai = N_ip1 - N_i - (c * (satellites[i][1] - satellites[i+1][1]))
        Xi = -(diff2[0]/N_ip1) + (diff[0]/N_i)
        Yi = -(diff2[1]/N_ip1) + (diff[1]/N_i)
        Zi = -(diff2[2]/N_ip1) + (diff[2]/N_i)

        s_i = satellites[i][2]
        s_ip1 = satellites[i+1][2]
        DXidx = DXiDx(N_ip1, s_ip1[0], N_i, s_i[0], x[0])
        DXidy = DYidx = DYiDx(N_ip1, s_ip1[1], s_ip1[0], N_i, s_i[1], s_i[0], x[0], x[1])
        DXidz = DZidx = DXiDz(N_ip1, s_ip1[0], s_ip1[2], N_i, s_i[0], s_i[2], x[0], x[2])
        DYidy = DYiDy(N_ip1, s_ip1[1], N_i, s_i[1], x[1])
        DYidz = DZidy = DYiDz(N_ip1, s_ip1[1], s_ip1[2], N_i, s_i[1], s_i[2], x[1], x[2])
        DZidz = DZiDz(N_ip1, s_ip1[2], N_i, s_i[2], x[2])

        DFdx += Ai * Xi
        DFdy += Ai * Yi
        DFdz += Ai * Zi

        J_inv += np.array([[(Xi * Xi) + (Ai * DXidx), (Xi * Yi) + (Ai * DXidy), (Xi * Zi) + (Ai * DXidz)],
                          [(Xi * Yi) + (Ai * DXidy), (Yi * Yi) + (Ai * DYidy), (Yi * Zi) + (Ai * DYidz)],
                          [(Xi * Zi) + (Ai * DXidz), (Yi * Zi) + (Ai * DYidz), (Zi * Zi) + (Ai * DZidz)]])

        pass

    F = np.array([2 * DFdx, 2 * DFdy, 2 * DFdz])
    J_inv = J_inv * 2
    return F, J_inv

# 
def calculateFirstOrderPartDeriv(x, satellites, startIndex, count) :
    vec = [mp.mpf(0.0),mp.mpf(0.0),mp.mpf(0.0)]
    for i in range(0,count-1) :
        index = startIndex + i
        n_i = magnitude(subVectors(satellites[index][2], x))
        n_iplus1 = magnitude(subVectors(satellites[index+1][2], x))
        a_i = n_iplus1 - n_i - c * (satellites[index][1] - satellites[index+1][1])
        df = addVectors(invScaleVector(-n_iplus1, subVectors(satellites[index+1], x)), invScaleVector(n_i, subVectors(satellites[index][2], x)))
        df = scaleVector(a_i, df)
        vec = addVectors(vec, df)

    vec = scaleVector(2, vec)
    return np.array([vec[0], vec[1], vec[2]])

#
def calculateSecOrderPartDeriv(x, satellites, startIndex, count) :
    for i in range(0,count-1) :
        index = startIndex + i
        n_i = magnitude(subVectors(satellites[index][2], x))
        n_iplus1 = magnitude(subVectors(satellites[index+1][2], x))
        DXidx = DXiDx(n_iplus1, satellites[index+1][2][0], n_i, satellites[index][2][0], x[0])
        DXidy = DYidx = DYiDx(n_iplus1, satellites[index+1][2][1], satellites[index+1][2][0], n_i, satellites[index][2][1], satellites[index][2][0], x[0], x[1])
        DXidz = DZidx = DXiDz(n_iplus1, satellites[index+1][2][0], satellites[index+1][2][2], n_i, satellites[index][2][0], satellites[index][2][2], x[0], x[2])
        DYidy = DYiDy(n_iplus1, satellites[index+1][2][1], n_i, satellites[index][2][1], x[1])
    return np.array([DXidx, DXidy, DXidz],
                    [DYidx, DYidy, DYidz],
                    [DZidx, DZidy, DZidz])

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
#xSi = None
#xSiPlus1 = None
#ySi = None
#ySiPlus1 = None
#zSi = None
#zSiPlus1 = None
#tSi = None
#tSiPlus1 = None
#NiPlus1 = None

# we make some generalizations here to make the first order partial derivatives easier to write
#Ni = magnitude(xSi - x)
#Ai = NiPlus1 - Ni - c(tSi - tSiPlus1)
#Xi = -(xSiPlus1 - x)/(NiPlus1) + (xSi - x)/Ni
#Yi = -(ySiPlus1 - y)/(NiPlus1) + (ySi - y)/Ni
#Zi = -(zSiPlus1 - z)/(NiPlus1) + (zSi - z)/Ni

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
    return -((ySiPlus1 - y) * (zSiPlus1 - z))/(NiPlus1^3) + ((ySi - y) * (zSi - z))/(Ni^3) #equivalent to DZiDy

def DZiDz(NiPlus1, zSiPlus1, Ni, zSi, z):
    return ((NiPlus1^2) - (zSiPlus1 - z)^2)/(NiPlus1^3) - ((Ni^2) - (zSi - z)^2)/(Ni^3)
    
main()
