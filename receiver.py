###########################################################
# 
# This program reads satellites locations/times and generates an approximation
# of a vehicle location on the globe.
# The mathematical approach came from homework 1. 
#
# Author(s): David Jones, Preston Malen, Leela Feaster
#
# This program was written for MATH5610 taught by professor Peter Alfield.
###########################################################

import sys
import mpmath as mp
import numpy as np
import os
import re

def main() :
    # create and set global variables
    global pi, s, c, r
    pi, s, c, r = getData()
    satellites, startIndicies = readSatelliteData()
    numPositions = len(startIndicies) - 1
    max_iterations = 200

    output_file = open(os.path.join(sys.path[0], 'receiver.log'), "w")

    # x_0, position in salt lake city at the first satellites signal time
    x_i = rotate(np.array([-1795225.28, -4477174.36, 4158593.45]), satellites[0][1])

    for i in range(0, numPositions) :

        startIndex = startIndicies[i]
        count = startIndicies[i+1] - startIndex

        # Newton's method to find the approx. vehicle position/time
        for k in range (0, max_iterations) :
            F, J_inv = calcF_Jinv(x_i, satellites, startIndex, count)
            try:
                s_i = np.linalg.solve(J_inv, -F)
            except np.linalg.LinAlgError as err:
                # matrix could be singular, just exit out early
                break
            x_i = x_i + s_i
            if(np.linalg.norm(s_i) < 0.01) :
                break
            pass

        # get t_v
        t_v = (np.linalg.norm(x_i - satellites[startIndex][2]) / c) + satellites[startIndex][1]

        # get geodesic coords
        h, lambda_d, lambda_m, lambda_s, phi_d, phi_m, phi_s, NS, EW = cartesianToGeodesic(rotate(x_i, -t_v), 2)
        sys.stdout.write("{:.2f} {} {} {:.3f} {} {} {} {:.2f} {} {:.2f}\n".format(t_v, phi_d, phi_m, phi_s, NS, lambda_d, lambda_m, lambda_s, EW, h))
        output_file.write("{:.2f} {} {} {:.3f} {} {} {} {:.2f} {} {:.2f}\n".format(t_v, phi_d, phi_m, phi_s, NS, lambda_d, lambda_m, lambda_s, EW, h))
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

    # only look at satellites that are above the horizon line
    validSatellites = list()
    x_dot = np.dot(x,x)
    for i in range(0, count) :
        s = satellites[startIndex + i]
        if(np.dot(x, s[2]) > x_dot) :
            validSatellites.append(s)
    validCount = len(validSatellites)

    # initialize local variables
    F = np.array([0.0, 0.0, 0.0])
    J_inv = np.array([[0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, 0.0]])

    # calculate F and inverse of the Jacobian at location x
    for i in range(0, validCount-1) :
        diff = validSatellites[i][2] - x
        diff2 = validSatellites[i+1][2] - x

        N_i = np.linalg.norm(diff)
        N_ip1 = np.linalg.norm(diff2)

        Ai = N_ip1 - N_i - (c * (validSatellites[i][1] - validSatellites[i+1][1]))
        Xi = -(diff2[0]/N_ip1) + (diff[0]/N_i)
        Yi = -(diff2[1]/N_ip1) + (diff[1]/N_i)
        Zi = -(diff2[2]/N_ip1) + (diff[2]/N_i)

        s_i = validSatellites[i][2]
        s_ip1 = validSatellites[i+1][2]
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
# assumes time is 0
# the seconds degree is rounded to the given digit
def cartesianToGeodesic(x, digit) :

    #calculate the height
    h = np.float64(np.linalg.norm(x) - r)

    # calculate psi
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
    elif(x[0] > 0.0 and x[1] < 0.0) :
        lamb = (2 * pi) + np.arctan2(x[1], x[0])
    elif(x[0] < 0.0) :
        lamb = pi + np.arctan2(x[1], x[0])

    if(lamb > pi) : 
        lamb = lamb - pi      

    # convert to degrees
    lamb = 180.0 * lamb / pi
    psi = 180.0 * psi / pi

    # get geodesic coords for lambda
    EW = int(np.sign(np.dot(np.array([0.0, 1.0, 0.0]), np.array([x[0], x[1], 0.0]))))
    lamb = 180.0 - lamb if (EW < 0) else lamb
    lambda_d = int(lamb)
    lamb = (lamb - lambda_d) * 60.0
    lambda_m = int(lamb)
    lamb -= lambda_m
    lambda_s = round(lamb * 60.0, digit)

    # account for rounding the seconds degree in lambda
    if abs(lambda_s - 60.0) < 0.00001 :
        lambda_s = 0.0
        lambda_m = lambda_m + 1

    if lambda_m - 60 == 0 :
        lambda_m = 0
        lambda_d = lambda_d + 1

    # get geodesic coords for psi
    NS = 1 if (psi >= 0) else -1
    psi = abs(psi)
    psi_d = int(psi)
    psi = (psi - psi_d) * 60.0
    psi_m = int(psi)
    psi = psi - psi_m
    psi_s = round(psi * 60.0, digit)

    # account for rounding the seconds degree in psi
    if abs(psi_s - 60.0) < 0.00001 :
        psi_s = 0.0
        psi_m = psi_m + 1
    
    if psi_m - 60 == 0 :
        psi_m = 0
        psi_d = psi_d + 1

    return h, lambda_d, lambda_m, lambda_s, psi_d, psi_m, psi_s, NS, EW

# rotates a vector along the z-axis
def rotate(u, t_v) :
    time_offset = (2*pi*t_v) / s
    v = np.array([(np.cos(time_offset) * u[0]) + (-np.sin(time_offset) * u[1]),
                  (np.sin(time_offset) * u[0]) + (np.cos(time_offset) * u[1]),
                   u[2]])
    return v

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
