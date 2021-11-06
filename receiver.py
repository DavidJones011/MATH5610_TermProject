import sys
import mpmath as mp
import numpy as np
import os
import re

def main() :
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
        x_0 = rotate(np.array([-1795225.28, -4477174.36, 4158593.45]), sattellites[startIndex][1])
        x_i = x_0

        # newtons method to find the vehicle location
        for k in range (0, max_iterations) :
            F, J_inv = calcF_Jinv(x_i, sattellites, startIndex, count)
            try:
                s_i = np.linalg.solve(J_inv, -F)
            except np.linalg.LinAlgError as err:
                break
            x_i = x_i + s_i
            if(abs(np.linalg.norm(s_i)) < 0.01) :
                break
            pass

        # get t_v
        t_v = (np.linalg.norm(x_i - sattellites[startIndex][2]) / c) + sattellites[startIndex][1]

        # get geodesic coords
        h, lambda_d, lambda_m, lambda_s, phi_d, phi_m, phi_s, NS, EW = cartesianToGeodesic(rotate(x_i, -t_v))
        sys.stdout.write("{} {} {} {} {} {} {} {} {} {}\n".format(t_v, phi_d, phi_m, phi_s, NS, lambda_d, lambda_m, lambda_s, EW, h))
        output_file.write("{} {} {} {} {} {} {} {} {} {}\n".format(t_v, phi_d, phi_m, phi_s, NS, lambda_d, lambda_m, lambda_s, EW, h))
        pass

    output_file.close()
    pass

# grab data from data.dat in local folder
def getData() :

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
        satellites.append(cur_satellite)

        # found a new set of satellites for a new position
        # store the starting index
        if cur_satellite[0] <= prevIndex : 
            startIndices.append(index)

        prevIndex = cur_satellite[0]
        index = index + 1

    startIndices.append(len(satellites))
    return satellites, startIndices

# grabs F(x) and the inverse Jacobian(x)
def calcF_Jinv(x, satellites, startIndex, count) :
    #initialize local variables
    F = np.array([0.0, 0.0, 0.0])
    J_inv = np.array([[0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0]])

    for i in range(0, count-1) :
        index = startIndex + i
        diff = satellites[index][2] - x
        diff2 = satellites[index+1][2] - x

        N_i = np.linalg.norm(diff)
        N_ip1 = np.linalg.norm(diff2)

        Ai = N_ip1 - N_i - (c * (satellites[index][1] - satellites[index+1][1]))
        Xi = -(diff2[0]/N_ip1) + (diff[0]/N_i)
        Yi = -(diff2[1]/N_ip1) + (diff[1]/N_i)
        Zi = -(diff2[2]/N_ip1) + (diff[2]/N_i)

        s_i = satellites[index][2]
        s_ip1 = satellites[index+1][2]
        DXidx = DXiDx(N_ip1, s_ip1[0], N_i, s_i[0], x[0])
        DXidy = DYiDx(N_ip1, s_ip1[1], s_ip1[0], N_i, s_i[1], s_i[0], x[0], x[1])
        DXidz = DXiDz(N_ip1, s_ip1[0], s_ip1[2], N_i, s_i[0], s_i[2], x[0], x[2])
        DYidy = DYiDy(N_ip1, s_ip1[1], N_i, s_i[1], x[1])
        DYidz = DYiDz(N_ip1, s_ip1[1], s_ip1[2], N_i, s_i[1], s_i[2], x[1], x[2])
        DZidz = DZiDz(N_ip1, s_ip1[2], N_i, s_i[2], x[2])

        F += np.array([Ai * Xi, Ai * Yi, Ai * Zi])

        J_inv += np.array([[(Xi * Xi) + (Ai * DXidx), (Xi * Yi) + (Ai * DXidy), (Xi * Zi) + (Ai * DXidz)],
                          [(Xi * Yi) + (Ai * DXidy), (Yi * Yi) + (Ai * DYidy), (Yi * Zi) + (Ai * DYidz)],
                          [(Xi * Zi) + (Ai * DXidz), (Yi * Zi) + (Ai * DYidz), (Zi * Zi) + (Ai * DZidz)]])

        pass

    F = F * 2.0
    J_inv = J_inv * 2.0
    return F, J_inv

# converts cartesian coordinate to geodesic
# assumes that time is 0
def cartesianToGeodesic(x) :
    #calculate the height
    h = np.float64(np.linalg.norm(x) - r)

    # calculate phi
    sq_dst = np.power(x[0], 2) + np.power(x[1], 2)
    psi = np.float64(0.0)
    if(sq_dst != 0.0) :
        psi = np.arctan2(x[2], np.sqrt(sq_dst))
    elif x[0] == x[1] and x[0] == 0.0:
        if(x[2] > 0.0) :
            psi = pi / 2.0
        elif(x[2] < 0.0) :
            psi = -pi / 2.0

    # calculate lambda
    lamb = np.float64(0.0)
    if(x[0] > 0.0 and x[1] > 0.0) :
        lamb = np.arctan2(x[1], x[0])
    elif(x[0] < 0.0) :
        lamb = pi + np.arctan2(x[1], x[0])
    elif(x[0] > 0.0 and x[1] < 0.0) :
        lamb = (2 * pi) + np.arctan2(x[1], x[0])

    # convert to degrees
    lamb = 180.0 * lamb / pi
    psi = 180.0 * psi / pi

    # get geodesic coords for lambda
    EW = int(np.sign(np.dot(np.array([1.0, 0.0, 0.0]), np.array([x[0], x[1], 0.0]))))
    lamb = 180.0 - lamb if (EW < 0) else lamb
    lambda_d = int(lamb)
    lamb = (lamb - lambda_d) * 60.0
    lambda_m = int(lamb)
    lamb -= lambda_m
    lambda_s = lamb * 60.0

    # get geodesic coords for phi
    NS = 1 if (psi >= 0) else -1
    psi_d = int(psi)
    psi = (psi - psi_d) * 60.0
    psi_m = int(psi)
    psi = psi - psi_m
    psi_s = psi * 60.0

    return h, lambda_d, lambda_m, lambda_s, psi_d, psi_m, psi_s, NS, EW

# rotates a vector along the z-axis
def rotate(u, t_v) :
    time_offset = (2*pi*t_v) / s
    v = np.array([(np.cos(time_offset) * u[0]) + (-np.sin(time_offset) * u[1]),
                  (np.sin(time_offset) * u[0]) + (np.cos(time_offset) * u[1]),
                   u[2]])
    return v

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
    return ((np.power(NiPlus1, 2) - np.power(xSiPlus1 - x, 2)) / np.power(NiPlus1, 3)) - ((np.power(Ni,2) - np.power(xSi - x, 2)) / np.power(Ni, 3))

def DYiDx(NiPlus1, ySiPlus1, xSiPlus1, Ni, ySi, xSi, x, y):
    return -(((ySiPlus1 - y) * (xSiPlus1 - x)) / np.power(NiPlus1, 3)) + (((ySi - y) * (xSi - x)) / np.power(Ni, 3)) #equivalent to DXiDy (partial of Xi with respect to y)

def DXiDz(NiPlus1, xSiPlus1, zSiPlus1, Ni, xSi, zSi, x, z):
    return -(((xSiPlus1 - x) * (zSiPlus1 - z)) / np.power(NiPlus1, 3)) + (((xSi - x) * (zSi - z)) / np.power(Ni, 3)) #equivalent to DZiDx

def DYiDy(NiPlus1, ySiPlus1, Ni, ySi, y):
    return ((np.power(NiPlus1, 2) - np.power(ySiPlus1 - y, 2)) / np.power(NiPlus1, 3)) - ((np.power(Ni, 2) - np.power(ySi - y, 2)) / np.power(Ni, 3))

def DYiDz(NiPlus1, ySiPlus1, zSiPlus1, Ni, ySi, zSi, y, z):
    return -(((ySiPlus1 - y) * (zSiPlus1 - z)) / np.power(NiPlus1, 3)) + (((ySi - y) * (zSi - z)) / np.power(Ni, 3)) #equivalent to DZiDy

def DZiDz(NiPlus1, zSiPlus1, Ni, zSi, z):
    return ((np.power(NiPlus1, 2) - np.power(zSiPlus1 - z, 2)) / np.power(NiPlus1, 3)) - ((np.power(Ni, 2) - np.power(zSi - z, 2)) / np.power(Ni, 3))
    
main()
