import sys
import math
import os
import re
from decimal import *

from enum import Enum
class GeoCoordinate(Enum) :
    LONGITUDE = 0
    LATITUDE = 1

def main() :
    # create global variables
    global pi, s, c, r
    # grab data from the data.dat file
    pi, s, c, r, satellites = getData()

    output_file = open(os.path.join(sys.path[0], 'satellite.log'), "w")

    for line in sys.stdin :
        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        #store the tokenized data into the correct data types
        t_v      = Decimal(tokens[0])
        psi_d    = int(tokens[1])
        psi_m    = int(tokens[2])
        psi_s    = Decimal(tokens[3])
        ns_value = int(tokens[4])
        lambda_d = int(tokens[5])
        lambda_m = int(tokens[6])
        lambda_s = Decimal(tokens[7])
        ew_value = int(tokens[8])
        altitude = Decimal(tokens[9])

        # get the cartesian coordinates of the vehicle
        latitude = degreeToRad(psi_d, psi_m, psi_s, ns_value, GeoCoordinate.LATITUDE)
        longitude = degreeToRad(lambda_d, lambda_m, lambda_s, ew_value, GeoCoordinate.LONGITUDE)
        x_v = sphericalToCartesian(r + altitude, longitude, latitude)
        dot_x_v = dotProduct(x_v,x_v)

        for i in range(0,24) :
            # use newtons method to find satellite position and time
            x_s, t_s = satelliteTimeAndLocOnSend(satellites, i, x_v, t_v)

            if dotProduct(x_v, x_s) > dot_x_v :       
                output_file.write("{} {} {:e} {:e} {:e}\n".format(i, t_s, x_s[0], x_s[1], x_s[2]))
            

    output_file.close()
    pass

# grab data from data.dat in local folder
def getData() :
    # initialize satellite data
    satellites = list()
    for i in range(0, 24):
        satellites.append([[0,0,0],[0,0,0],0,0,0])

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
            pi = Decimal(tokens[0])
        elif(tokens[2] == 's') :
            s = Decimal(tokens[0])
        elif(tokens[2] == 'c') :
            c = Decimal(tokens[0])
        elif(tokens[2] == 'R') :
            r = Decimal(tokens[0])
        elif(tokens[2] == 'u1') :
            satellites[int(tokens[5])][0][0] = Decimal(tokens[0])
        elif(tokens[2] == 'u2') :
            satellites[int(tokens[5])][0][1] = Decimal(tokens[0])
        elif(tokens[2] == 'u3') :
            satellites[int(tokens[5])][0][2] = Decimal(tokens[0])
        elif(tokens[2] == 'v1') :
            satellites[int(tokens[5])][1][0] = Decimal(tokens[0])
        elif(tokens[2] == 'v2') :
            satellites[int(tokens[5])][1][1] = Decimal(tokens[0])
        elif(tokens[2] == 'v3') :
            satellites[int(tokens[5])][1][2] = Decimal(tokens[0])
        elif(tokens[2] == 'periodicity') :
            satellites[int(tokens[5])][2] = Decimal(tokens[0])
        elif(tokens[2] == 'altitude') :
            satellites[int(tokens[5])][3] = Decimal(tokens[0])
        elif(tokens[2] == 'phase') :
            satellites[int(tokens[5])][4] = Decimal(tokens[0])

    data_file.close()
    return pi, s, c, r, satellites

# converts degree (geographic) to radians
def degreeToRad(angle, angle_minute, angle_second, dir_value, geo_coordinate) :
    absolute_angle = Decimal(angle) + (Decimal(angle_minute) / 60) + (angle_second / 3600)

    # convert the angle to be within (0,360)
    if(geo_coordinate == GeoCoordinate.LONGITUDE) :
        absolute_angle =  ((1 - ((1 + dir_value) * Decimal(0.5))) * 360) + (dir_value * absolute_angle)
        
    # convert the angle to be within (0, 180)   
    elif(geo_coordinate == GeoCoordinate.LATITUDE) :
        absolute_angle = 90 - (dir_value * absolute_angle)
     
    return absolute_angle * (pi / 180)

# converts radians to degrees (geographic) #NOT IN USE
def radToDegree(angle, geo_coordinate) :
    absolute_angle = angle * 180.0 / pi
    dir_value = 1

    # convert the range from (0, 360) to (0, 180)
    # also gets the EW value
    if(geo_coordinate == GeoCoordinate.LONGITUDE) :
        dir_value = int(-1 if (absolute_angle > pi) else 1)
        absolute_angle = (absolute_angle - ((1 - ((1 + dir_value) * 0.5)) * 360.0)) * dir_value
    
    # converts the range from (0, 180) to (0,90)
    # also gets the NS value
    if(geo_coordinate == GeoCoordinate.LATITUDE) :
        dir_value = int(-1 if (absolute_angle > 90) else 1)
        absolute_angle = (absolute_angle - 90.0) * -dir_value

    # calculates the angle in degree, minutes, and seconds
    angle = int(math.floor(absolute_angle))
    absolute_angle = (absolute_angle - angle) * 60.0
    angle_min = int(math.floor(absolute_angle))
    angle_sec = float((absolute_angle - angle_min) * 60.0)

    return [angle, angle_min, angle_sec, dir_value]

# get cartesian coords from spherical coords
def sphericalToCartesian(radius, angle_theta, angle_phi) :
    x = radius * Decimal(math.cos(angle_theta) * math.sin(angle_phi))
    y = radius * Decimal(math.sin(angle_theta) * math.sin(angle_phi))
    z = radius * Decimal(math.cos(angle_phi))
    return [x, y, z]

# get geographic coords from cartesian coords #NOT IN USE
def cartesianToGeographic(x) :
    mag = magnitude(x)
    altitude = mag - r
    mag2 = Decimal.sqrt(x[0] * x[0] + x[1] * x[1])
    theta = math.atan(x[1]/ x[0]) + pi
    phi = (pi / 2) - math.atan(x[2] / mag2)
    return [theta, phi, altitude]

# get the cartesian coordinates of satellite
def getSatelliteLocation(satellites, index, t) :
    angle = (2 * pi * t) / satellites[index][2]
    coeff = r + satellites[index][3]
    u = scaleVector(math.cos(angle + satellites[index][4]), satellites[index][0])
    v = scaleVector(math.sin(angle + satellites[index][4]), satellites[index][1])
    return scaleVector(coeff, addVectors(u,v))

# get the satellites velocity
def getSatelliteVelocity(satellites, index, t) :
    coeff = (2 * pi * (r + satellites[index][3])) / satellites[index][2]
    angle = (2 * pi * t) / satellites[index][2]
    u = scaleVector(-1 * math.sin(angle + satellites[index][4]), satellites[index][0])
    v = scaleVector(math.cos(angle + satellites[index][4]), satellites[index][1])
    return scaleVector(coeff, addVectors(u,v))

# get the satellite location and time to send the signal at x_v at time t
def satelliteTimeAndLocOnSend(satellites, index, x_v, t_v) :
    best_t = Decimal(t_v)
    x0_s = getSatelliteLocation(satellites, index, t_v)
    prev_t = Decimal(0)
    cur_t = t_v - (magnitude(subVectors(x0_s, x_v)) / c)
    threshold = Decimal(0.01) / c

    # we want to have number of iterations in case NM diverges
    while Decimal.copy_abs(cur_t - prev_t) > threshold : 
        cur_x = getSatelliteLocation(satellites, index, cur_t)
        f = magnitude(subVectors(cur_x, x_v)) - (c * (t_v - cur_t))
        f_prime = (2 * dotProduct(subVectors(cur_x, x_v), getSatelliteVelocity(satellites, index, cur_t))) + (2 * (c * c) * (t_v - cur_t))

        # in the case the f' is not valid
        if(math.isinf(f_prime)) :
            break

        prev_t = cur_t
        cur_t = prev_t - (f / f_prime)
        best_t = t_v + (cur_t - t_v)

    return getSatelliteLocation(satellites, index, best_t), best_t

# returns a vector that is scaled
def scaleVector(scale, vec) :
    x = Decimal(vec[0] * Decimal(scale))
    y = Decimal(vec[1] * Decimal(scale))
    z = Decimal(vec[2] * Decimal(scale))
    return [x, y, z]

# returns the magnitude of u
def magnitude(u):
    return Decimal.sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]))

# returns the dot product of two vectors u and v
def dotProduct(u, v) :
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2])

# addes two vectors
def addVectors(u,v) :
    return [u[0] + v[0], u[1] + v[1], u[2] + v[2]]

def subVectors(u,v) :
    return [u[0] - v[0], u[1] - v[1], u[2] - v[2]]

main()