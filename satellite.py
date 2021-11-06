###########################################################
# 
# Description
#
# Author(s): David Jones, Preston , Leela
#
###########################################################

import sys
import os
import re
import numpy as np

def main() :
    # create and set global variables
    global pi, s, c, r
    pi, s, c, r, satellites = getData()
    vehicles = getVehcileData()

    output_file = open(os.path.join(sys.path[0], 'satellite.log'), "w")

    for v in vehicles :
        # get the cartesian coordinates of the vehicle
        latitude = degreeToRad(v[1], v[2], v[3], v[4])
        longitude = degreeToRad(v[5], v[6], v[7], v[8])
        x_v = sphericalToCartesian(r + v[9], longitude, latitude)
        x_v = rotate(x_v, v[0])

        dot_x_v = np.dot(x_v,x_v)

        for i in range(0,24) :
            # use newtons method to find satellite position and time
            x_s, t_s = satelliteTimeAndLocOnSend(satellites, i, x_v, v[0])

            if np.dot(x_v, x_s) > dot_x_v :
                output_file.write("{} {} {} {} {}\n".format(i, t_s, x_s[0], x_s[1], x_s[2]))
                sys.stdout.write("{} {} {} {} {}\n".format(i, t_s, x_s[0], x_s[1], x_s[2]))
            

    output_file.close()
    pass

# grab the vehicle position/time data
def getVehcileData() :

    vehicles = list()

    for line in sys.stdin :
        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        #store the tokenized data into the correct data types
        cur_vehicle = [0.0, 0, 0, 0.0, 0, 0, 0, 0.0, 0, 0.0]
        cur_vehicle[0] = np.float64(tokens[0])    # t_v
        cur_vehicle[1] = int(tokens[1])           # psi_d
        cur_vehicle[2] = int(tokens[2])           # psi_m
        cur_vehicle[3] = np.float64(tokens[3])    # psi_s
        cur_vehicle[4] = int(tokens[4])           # NS
        cur_vehicle[5] = int(tokens[5])           # lambda_d
        cur_vehicle[6] = int(tokens[6])           # lambda_m
        cur_vehicle[7] = np.float64(tokens[7])    # lambda_s
        cur_vehicle[8] = int(tokens[8])           # EW
        cur_vehicle[9] = np.float64(tokens[9])    # altitude
        vehicles.append(cur_vehicle)

    return vehicles

# grab data from data.dat in local folder
def getData() :
    # initialize satellite data
    satellites = list()
    for i in range(0, 24):
        satellites.append([np.array([0.0,0.0,0.0]), np.array([0.0,0.0,0.0]), 0.0, 0.0, 0.0])

    data_file = open(os.path.join(sys.path[0], 'data.dat'), "r")
    for line in data_file.readlines() :

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
        elif(tokens[2] == 'u1') :
            satellites[int(tokens[5])][0][0] = np.float64(tokens[0])
        elif(tokens[2] == 'u2') :
            satellites[int(tokens[5])][0][1] = np.float64(tokens[0])
        elif(tokens[2] == 'u3') :
            satellites[int(tokens[5])][0][2] = np.float64(tokens[0])
        elif(tokens[2] == 'v1') :
            satellites[int(tokens[5])][1][0] = np.float64(tokens[0])
        elif(tokens[2] == 'v2') :
            satellites[int(tokens[5])][1][1] = np.float64(tokens[0])
        elif(tokens[2] == 'v3') :
            satellites[int(tokens[5])][1][2] = np.float64(tokens[0])
        elif(tokens[2] == 'periodicity') :
            satellites[int(tokens[5])][2] = np.float64(tokens[0])
        elif(tokens[2] == 'altitude') :
            satellites[int(tokens[5])][3] = np.float64(tokens[0])
        elif(tokens[2] == 'phase') :
            satellites[int(tokens[5])][4] = np.float64(tokens[0])

    data_file.close()
    return pi, s, c, r, satellites

# converts degree (geodesic) to radians
def degreeToRad(angle, angle_minute, angle_second, dir_value) :
    absolute_angle = angle + (np.float64(angle_minute) / 60.0) + (np.float64(angle_second) / 3600.0)
    return absolute_angle * (pi / 180.0) * dir_value

# get cartesian coords from spherical coords
def sphericalToCartesian(radius, angle_lambda, angle_phi) :
    x = radius * np.cos(angle_lambda) * np.cos(angle_phi)
    y = radius * np.sin(angle_lambda) * np.cos(angle_phi)
    z = radius * np.sin(angle_phi)
    return np.array([x, y, z])

# get the cartesian coordinates of satellite
def getSatelliteLocation(satellites, index, t) :
    angle = 2.0 * pi * t / satellites[index][2]
    coeff = r + satellites[index][3]
    u = satellites[index][0] * np.cos(angle + satellites[index][4])
    v = satellites[index][1] * np.sin(angle + satellites[index][4])
    return (u + v) * coeff

# get the satellites velocity
def getSatelliteVelocity(satellites, index, t) :
    coeff = (2.0 * pi * (r + satellites[index][3])) / satellites[index][2]
    angle = (2.0 * pi * t) / satellites[index][2]
    u = satellites[index][0] * -np.sin(angle + satellites[index][4])
    v = satellites[index][1] * np.cos(angle + satellites[index][4])
    return (u + v) * coeff

# get the satellite location and time to send the signal at x_v at time t
def satelliteTimeAndLocOnSend(satellites, index, x_v, t_v) :
    x0_s = getSatelliteLocation(satellites, index, t_v)
    # calculate t0
    prev_t = t_v - (np.linalg.norm(x0_s - x_v) / c)
    threshold = 0.01 / c
    max_iterations = 200
    best_t = 0

    # we want to have number of iterations in case NM diverges
    for i in range(0,max_iterations) :
        pos = getSatelliteLocation(satellites, index, prev_t)

        # get f(x)
        diff = pos - x_v
        f = np.dot(diff, diff) - (np.power(c, 2) * (np.power(t_v - prev_t, 2)))
        # get f'(x)
        f_prime = (2 * np.dot(diff, getSatelliteVelocity(satellites, index, prev_t))) + (2.0 * np.power(c, 2) * (t_v - prev_t))
        cur_t = prev_t - (f / f_prime)

        if np.abs(cur_t - prev_t) < threshold :
            best_t = cur_t
            break

        prev_t = cur_t

    return getSatelliteLocation(satellites, index, best_t), best_t

# returns a vector that is scaled
#def scaleVector(scale, vec) :
    x = vec[0] * scale
    y = vec[1] * scale
    z = vec[2] * scale
    return [x, y, z]

# returns the magnitude of u
#def magnitude(u):
    return mp.sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]))

# returns the dot product of two vectors u and v
#def dotProduct(u, v) :
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2])

# adds two vectors
#def addVectors(u,v) :
    return [u[0] + v[0], u[1] + v[1], u[2] + v[2]]

# subtracts two vectors
#def subVectors(u,v) :
    return [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

# rotates a point along the z-axis
def rotate(u, t_v) :
    time_offset = (2*pi*t_v) / s
    v = np.array([(np.cos(time_offset) * u[0]) + (-np.sin(time_offset) * u[1]),
         (np.sin(time_offset) * u[0]) + (np.cos(time_offset) * u[1]),
         (u[2])])
    return v

main()