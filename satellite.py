import sys
import math
import os
import re

from enum import Enum
class GeoCoordinate(Enum) :
    LONGITUDE = 0
    LATITUDE = 1

from enum import Enum
class Axis(Enum) :
    X = 0
    Y = 1
    Z = 2

SPEED_OF_LIGHT = float(2.997924580000000000E+08)
EARTH_RADIUS = float(6.367444500000000000E+06)
SIDEREAL_DAY_SEC = float(8.616408999999999651E+04)
PI = float(3.141592653589793116E+00)

def main() :
    # open the data file to retrieve info
    data_file = open(os.path.join(sys.path[0], 'data.dat'), "r")
    for line in data_file.readlines() :
        if line == '':
            break
        
        tokens = re.split('\s|,', line)

        

        if(len(tokens) == 0) :
            break

        
        
    data_file.close()

    for line in sys.stdin :
        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        #store the tokenized data into the correct data types
        t_v  = float(tokens[0])
        psi_d    = int(tokens[1])
        psi_m    = int(tokens[2])
        psi_s    = float(tokens[3])
        ns_value = int(tokens[4])
        lambda_d = int(tokens[5])
        lambda_m = int(tokens[6])
        lambda_s = float(tokens[7])
        ew_value = int(tokens[8])
        altitude   = float(tokens[9])

        #lat_rad = degreeToRad(psi_d, psi_m, psi_s, ns_value, GeoCoordinate.LATITUDE)
        #lon_rad = degreeToRad(lambda_d, lambda_m, lambda_s, ew_value, GeoCoordinate.LONGITUDE)
        #x_v = getCartesianCoords(EARTH_RADIUS + altitude, lon_rad, lat_rad)

    newvector = rotateVectorAroundAxis([5,0,0], math.pi/2.0, Axis.Y)
    newvector = normalize(newvector)

    # writes the info into a file
    output_file = open(os.path.join(sys.path[0], 'satellite.log'), "w")
    output_file.write("<{},{},{}>\n".format(newvector[0], newvector[1], newvector[2]))
    output_file.close()
    pass

def degreeToRad(angle, angle_minute, angle_second, dir_value, geo_coordinate) :
    absolute_angle = angle + (angle_minute / 60.0) + (angle_second / 3600.0)

    # convert the angle to be within (0,360)
    if(geo_coordinate == GeoCoordinate.LONGITUDE) :
        absolute_angle =  ((1 - ((1 + dir_value) * 0.5)) * 360.0) + (dir_value * absolute_angle) 
        
    # convert the angle to be within (0, 180)   
    elif(geo_coordinate == GeoCoordinate.LATITUDE) :
        absolute_angle = 90.0 - (dir_value * absolute_angle)
     
    return absolute_angle * PI / 180.0

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

def getCartesianCoords(radius, angle_theta, angle_phi) :
    x = radius * math.cos(angle_theta)*math.sin(angle_phi)
    y = radius * math.sin(angle_theta)*math.sin(angle_phi)
    z = radius * math.cos(angle_phi)
    return [x, y, z]

def normalize(u) :
    mag = magnitude(u)
    return [u[0] / mag, u[1] / mag, u[2] / mag]

def magnitude(u):
    return math.sqrt((u[0] * u[0]) + (u[1] * u[1]) + (u[2] * u[2]))

def dotProduct(u, v) :
    return (u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2])

def rotateVectorAroundAxis(u, angle, axis) :
    v = u
    if axis == Axis.X :
        v = [u[0] 
            ,u[1] * math.cos(angle) - u[2] * math.sin(angle)
            ,u[1] * math.sin(angle) - u[2] * math.cos(angle)]
    if axis == Axis.Y :
        v = [u[0] * math.cos(angle) + u[2] * math.cos(angle)
            ,u[1]
            ,-u[0] * math.sin(angle) + u[2] * math.cos(angle)]
    if axis == Axis.Z :
        v = [u[0] * math.cos(angle) - u[1] * math.sin(angle)
            ,u[0] * math.sin(angle) + u[1] * math.cos(angle)
            ,u[2]]
    return v

main()