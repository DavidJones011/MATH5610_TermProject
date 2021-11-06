import sys
import os
import re
import numpy as np

def main() :

    line = sys.stdin.readline()
    tokens = re.split('\s|,', line)

    #remove whitespace
    while("" in tokens) :
        tokens.remove("")

    if(len(tokens) != 2) :
        sys.stdout.write("Input was invalid.\n")
        return

    thisVehicles = getVehcileData(tokens[0])
    otherVehicles = getVehcileData(tokens[1])

    if(len(thisVehicles) != len(otherVehicles)) :
        sys.stdout.write("Failed! Number of vehicle postions/times mismatched.\n")
        return


    success = True
    for i in range(0, len(thisVehicles)) :

        tVehicle = thisVehicles[i]
        oVehicle = otherVehicles[i]
        
        t_lamb = degreeToRad(tVehicle[5], tVehicle[6], tVehicle[7], tVehicle[8])
        t_psi = degreeToRad(tVehicle[1], tVehicle[2], tVehicle[3], tVehicle[4])
        o_lamb = degreeToRad(oVehicle[5], oVehicle[6], oVehicle[7], oVehicle[8])
        o_psi = degreeToRad(oVehicle[1], oVehicle[2], oVehicle[3], oVehicle[4]) 

        tLoc = sphericalToCartesian(tVehicle[9], t_lamb, t_psi)
        oLoc = sphericalToCartesian(tVehicle[9], o_lamb, o_psi)
        dist = round(np.linalg.norm(tLoc - oLoc), 2)

        if(dist > 0.001) :
            sys.stdout.write("Vehicle locations differ further than a centimeter! Line: {}\n".format(i))
            success = False
            break

        tTime = round(thisVehicles[i][0], 15)
        oTime = round(otherVehicles[i][0], 15)

        # could test for accuracy in seconds
        if(tTime - oTime > 0.001) :
            sys.stdout.write("Vehcile position times differ too much! Line: {}\n".format(i))
            success = False
            break

        pass

    if success :
        sys.stdout.write("Success! The results were close.\n")
    else :
        sys.stdout.write("Failed! :(\n")

    pass

# grab the vehicle position/time data
def getVehcileData(fileName) :
    data_file = open(os.path.join(sys.path[0], fileName), "r")
    vehicles = list()

    for line in data_file.readlines() :

        if line == '':
            break

        tokens = re.split('\s|,', line)

        if(len(tokens) == 0) :
            break

        # store the tokenized data into the correct data types
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

    data_file.close()
    return vehicles

# converts degrees to radians
def degreeToRad(angle, angle_minute, angle_second, dir_value) :
    absolute_angle = angle + (np.float64(angle_minute) / 60.0) + (np.float64(angle_second) / 3600.0)
    return absolute_angle * (np.pi / 180.0) * dir_value

# get cartesian coords from spherical coords
def sphericalToCartesian(radius, angle_lambda, angle_phi) :
    x = radius * np.cos(angle_lambda) * np.cos(angle_phi)
    y = radius * np.sin(angle_lambda) * np.cos(angle_phi)
    z = radius * np.sin(angle_phi)
    return np.array([x, y, z])

main()