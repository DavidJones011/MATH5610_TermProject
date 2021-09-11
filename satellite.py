import sys
import math
import os
import re
import decimal

from enum import Enum
class GeoCoordinate(Enum) :
    LONGITUDE = 0
    LATITUDE = 1

# global variables
pi = 0
s = 0 
c = 0
r = 0

def main() :
    pi, s, c, r, satellites = getData()

    for line in sys.stdin :
        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        #store the tokenized data into the correct data types
        t_v      = float(tokens[0])
        psi_d    = int(tokens[1])
        psi_m    = int(tokens[2])
        psi_s    = float(tokens[3])
        ns_value = int(tokens[4])
        lambda_d = int(tokens[5])
        lambda_m = int(tokens[6])
        lambda_s = float(tokens[7])
        ew_value = int(tokens[8])
        altitude = float(tokens[9])

        # get the cartesian coordinates of the vehicle
        latitude = degreeToRad(psi_d, psi_m, psi_s, ns_value, GeoCoordinate.LATITUDE)
        longitude = degreeToRad(lambda_d, lambda_m, lambda_s, ew_value, GeoCoordinate.LONGITUDE)
        x_v = sphericalToCartesian(r + altitude, longitude, latitude)

        # grab the satellites that are considered above the horizon plane
        sat_whitelist = list()
        for i in range(0,24) :
            x_s = getSatelliteCartesian(satellites, i, t_v)
            dot = dotProduct(normalize(x_v), normalize(x_s, -x_v))
            if(dot < 0) :
                continue
            sat_whitelist.append(i)

        for i in sat_whitelist :
            #
            # project satellite onto earth and find difference in time
             
        #output_file = open(os.path.join(sys.path[0], 'satellite.log'), "w")
        #output_file.write("<{},{},{}>\n".format(newvector[0], newvector[1], newvector[2]))
        #output_file.close()
    pass

# grab data from data.dat in local folder
def getData() :
    # initialize satellite data
    satellites = list()
    for i in range(0, 24):
        offset = decimal.Decimal((i % 4) * (math.pi / 2))
        satellites.append([[0,0,0],[0,0,0],0,0,0, offset])

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
            pi = decimal.Decimal(tokens[0])
        elif(tokens[2] == 's') :
            s = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'c') :
            c = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'R') :
            r = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'u1') :
            satellites[int(tokens[5])][0][0] = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'u2') :
            satellites[int(tokens[5])][0][1] = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'u3') :
            satellites[int(tokens[5])][0][2] = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'v1') :
            satellites[int(tokens[5])][1][0] = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'v2') :
            satellites[int(tokens[5])][1][1] = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'v3') :
            satellites[int(tokens[5])][1][2] = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'periodicity') :
            satellites[int(tokens[5])][2] = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'altitude') :
            satellites[int(tokens[5])][3] = decimal.Decimal(tokens[0])
        elif(tokens[2] == 'phase') :
            satellites[int(tokens[5])][4] = decimal.Decimal(tokens[0])

    data_file.close()
    return pi, s, c, r, satellites

# converts degree (geographic) to radians
def degreeToRad(angle, angle_minute, angle_second, dir_value, geo_coordinate) :
    absolute_angle = angle + (angle_minute / 60.0) + (angle_second / 3600.0)

    # convert the angle to be within (0,360)
    if(geo_coordinate == GeoCoordinate.LONGITUDE) :
        absolute_angle =  ((1 - ((1 + dir_value) * 0.5)) * 360.0) + (dir_value * absolute_angle) 
        
    # convert the angle to be within (0, 180)   
    elif(geo_coordinate == GeoCoordinate.LATITUDE) :
        absolute_angle = 90.0 - (dir_value * absolute_angle)
     
    return absolute_angle * PI / 180.0

# converts radians to degrees (geographic)
def radToDegree(angle, geo_coordinate) :
    absolute_angle = angle * 180.0 / PI
    dir_value = 1

    # convert the range from (0, 360) to (0, 180)
    # also gets the EW value
    if(geo_coordinate == GeoCoordinate.LONGITUDE) :
        dir_value = int(-1 if (absolute_angle > PI) else 1)
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
    x = radius * math.cos(angle_theta)*math.sin(angle_phi)
    y = radius * math.sin(angle_theta)*math.sin(angle_phi)
    z = radius * math.cos(angle_phi)
    return [x, y, z]

# get geographic coords from cartesian coords
def cartesianToGeographic(x) :
    mag = magnitude(x)
    altitude = mag - r
    mag2 = math.sqrt(x[0] * x[0] + x[1] * x[1])
    theta = math.atan(x[1]/ x[0]) + math.pi
    phi = (math.pi / 2) - math.atan(x[2] / mag2)
    return [theta, phi, altitude]

# get the cartesian coordinates of satellite
def getSatelliteCartesian(satellites, index, t) :
    # out of bounds
    if(index >= satellites.length and index < 0) :
        return [0, 0, 0]   
    angle = (2 * pi * t) / satellites[2]
    u = scaleVector((r + satellites[index][3]) * math.cos(angle + satellites[index][5]), satellites[index][0])
    v = scaleVector((r + satellites[index][3]) * math.sin(angle + satellites[index][5]), satellites[index][1])
    return addVectors(u,v)

# returns a vector that is scaled
def scaleVector(scale, vec) :
    x = vec * scale
    y = vec * scale
    z = vec * scale
    return [x, y, z]

# returns a normalized vector of u
def normalize(u) :
    return scaleVector(1/magnitude(u), u)

# returns the magnitude of u
def magnitude(u):
    return math.sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]))

# returns the dot product of two vectors u and v
def dotProduct(u, v) :
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2])

# addes two vectors
def addVectors(u,v) :
    return [u[0] + v[0], u[1] + v[1], u[2] + v[2]]

#rotation in xyz order
def rotateVectorAroundAxis(vec, roll, pitch, yaw) :
    x = vec[0] * (math.cos(yaw) * math.cos(pitch)) + vec[1] * ((math.cos(yaw) * math.sin(pitch) * math.sin(roll)) - (math.sin(yaw) * math.cos(roll))) + vec[2] * ((math.cos(yaw) * math.sin(pitch) * math.cos(roll)) + (math.sin(yaw) * math.sin(roll)))
    y = vec[0] * (math.sin(yaw) * math.cos(pitch)) + vec[1] * ((math.sin(yaw) * math.sin(pitch) * math.sin(roll)) + (math.cos(yaw) * math.cos(roll))) + vec[2] * ((math.sin(yaw) * math.sin(pitch) * math.cos(roll)) - (math.cos(yaw) * math.sin(roll)))
    z = -vec[0] * math.sin(pitch) + vec[1] * math.cos(pitch) * math.sin(roll) + vec[2] * math.cos(pitch) * math.cos(roll)
    return [x,y,z]

# call the program
main()