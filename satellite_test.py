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

    thisSatellites = getSatelliteData(tokens[0])
    otherSatellites = getSatelliteData(tokens[1])

    if(len(thisSatellites) != len(otherSatellites)) :
        sys.stdout.write("Failed! Number of satellites mismatched.\n")
        return


    success = True
    for i in range(0, len(thisSatellites)) :

        tIndex = thisSatellites[i][0]
        oIndex = otherSatellites[i][0]

        if(tIndex != oIndex) :
            sys.stdout.write("Satellite index mismatch! {} and {}\n".format(tIndex, oIndex))
            success = False
            break

        tLoc = thisSatellites[i][2]
        oLoc = otherSatellites[i][2]

        if(round(np.linalg.norm(tLoc - oLoc), 2) > 0.001) :
            sys.stdout.write("Satellite locations differ further than a centimeter! {} and {}\n".format(tIndex, oIndex))
            success = False
            break

        tTime = round(thisSatellites[i][1], 2)
        oTime = round(otherSatellites[i][1], 2)

        # could test for accuracy in seconds
        if(tTime - oTime > 0.001) :
            sys.stdout.write("Satellite times differ too much! {} and {}\n".format(tIndex, oIndex))
            success = False
            break

        pass

    if success :
        sys.stdout.write("Success! The results were close.\n")
    else :
        sys.stdout.write("Failed! :(\n")

    pass

def getSatelliteData(fileName) :
    data_file = open(os.path.join(sys.path[0], fileName), "r")
    satellites = list()
    i = 0
    for line in data_file.readlines() :

        if line == '':
            break

        tokens = re.split('\s|,', line)

        # remove empty strings from the list of tokens
        while("" in tokens) :
            tokens.remove("")

        if(len(tokens) == 0) :
            break

        s = [0, 0, [0,0,0]]
        s[0]    = int(tokens[0])
        s[1]    = np.float64(tokens[1])
        s[2] = np.array([np.float64(tokens[2]), np.float64(tokens[3]), np.float64(tokens[4])])
        satellites.append(s)
        i = i + 1

    data_file.close()

    return satellites

main()