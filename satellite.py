import sys
import os
import re
import mpmath as mp

def main() :
    # set the decimal place for calculations
    mp.mp.dps = 17
    # create and set global variables
    global pi, s, c, r
    pi, s, c, r, satellites = getData()

    output_file = open(os.path.join(sys.path[0], 'satellite.log'), "w")

    for line in sys.stdin :
        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        #store the tokenized data into the correct data types
        t_v      = mp.mpf(tokens[0])
        psi_d    = int(tokens[1])
        psi_m    = int(tokens[2])
        psi_s    = mp.mpf(tokens[3])
        ns_value = int(tokens[4])
        lambda_d = int(tokens[5])
        lambda_m = int(tokens[6])
        lambda_s = mp.mpf(tokens[7])
        ew_value = int(tokens[8])
        altitude = mp.mpf(tokens[9])

        # get the cartesian coordinates of the vehicle
        latitude = degreeToRad(psi_d, psi_m, psi_s, ns_value)
        longitude = degreeToRad(lambda_d, lambda_m, lambda_s, ew_value)
        x_v = sphericalToCartesian(r + altitude, longitude, latitude)
        x_v = rotate(x_v, t_v)

        dot_x_v = dotProduct(x_v,x_v)

        for i in range(0,24) :
            # use newtons method to find satellite position and time
            x_s, t_s = satelliteTimeAndLocOnSend(satellites, i, x_v, t_v)

            if dotProduct(x_v, x_s) > dot_x_v :
                output_file.write("{} {} {} {} {}\n".format(i, t_s, x_s[0], x_s[1], x_s[2]))
                sys.stdout.write("{} {} {} {} {}\n".format(i, t_s, x_s[0], x_s[1], x_s[2]))
            

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
            pi = mp.mpf(tokens[0])
        elif(tokens[2] == 's') :
            s = mp.mpf(tokens[0])
        elif(tokens[2] == 'c') :
            c = mp.mpf(tokens[0])
        elif(tokens[2] == 'R') :
            r = mp.mpf(tokens[0])
        elif(tokens[2] == 'u1') :
            satellites[int(tokens[5])][0][0] = mp.mpf(tokens[0])
        elif(tokens[2] == 'u2') :
            satellites[int(tokens[5])][0][1] = mp.mpf(tokens[0])
        elif(tokens[2] == 'u3') :
            satellites[int(tokens[5])][0][2] = mp.mpf(tokens[0])
        elif(tokens[2] == 'v1') :
            satellites[int(tokens[5])][1][0] = mp.mpf(tokens[0])
        elif(tokens[2] == 'v2') :
            satellites[int(tokens[5])][1][1] = mp.mpf(tokens[0])
        elif(tokens[2] == 'v3') :
            satellites[int(tokens[5])][1][2] = mp.mpf(tokens[0])
        elif(tokens[2] == 'periodicity') :
            satellites[int(tokens[5])][2] = mp.mpf(tokens[0])
        elif(tokens[2] == 'altitude') :
            satellites[int(tokens[5])][3] = mp.mpf(tokens[0])
        elif(tokens[2] == 'phase') :
            satellites[int(tokens[5])][4] = mp.mpf(tokens[0])

    data_file.close()
    return pi, s, c, r, satellites

# converts degree (geodesic) to radians
def degreeToRad(angle, angle_minute, angle_second, dir_value) :
    absolute_angle = angle + mp.mpf(angle_minute / 60.0) + mp.mpf(angle_second / 3600.0)
    return absolute_angle * (pi / 180.0) * dir_value

# get cartesian coords from spherical coords
def sphericalToCartesian(radius, angle_lambda, angle_phi) :
    x = radius * mp.cos(angle_lambda) * mp.cos(angle_phi)
    y = radius * mp.sin(angle_lambda) * mp.cos(angle_phi)
    z = radius * mp.sin(angle_phi)
    return [x, y, z]

# get the cartesian coordinates of satellite
def getSatelliteLocation(satellites, index, t) :
    angle = mp.mpf((2.0 * pi * t) / satellites[index][2])
    coeff = r + satellites[index][3]
    u = scaleVector(mp.cos(angle + satellites[index][4]), satellites[index][0])
    v = scaleVector(mp.sin(angle + satellites[index][4]), satellites[index][1])
    return scaleVector(coeff, addVectors(u,v))

# get the satellites velocity
def getSatelliteVelocity(satellites, index, t) :
    coeff = (2.0 * pi * (r + satellites[index][3])) / satellites[index][2]
    angle = (2.0 * pi * t) / satellites[index][2]
    u = scaleVector(-1.0 * mp.sin(angle + satellites[index][4]), satellites[index][0])
    v = scaleVector(mp.cos(angle + satellites[index][4]), satellites[index][1])
    return scaleVector(coeff, addVectors(u,v))

# get the satellite location and time to send the signal at x_v at time t
def satelliteTimeAndLocOnSend(satellites, index, x_v, t_v) :
    x0_s = getSatelliteLocation(satellites, index, t_v)
    # calculate t0
    prev_t = t_v - (magnitude(subVectors(x0_s, x_v)) / c)
    threshold = 0.01 / c
    max_iterations = 200
    best_t = 0

    # we want to have number of iterations in case NM diverges
    for i in range(0,max_iterations) :
        pos = getSatelliteLocation(satellites, index, prev_t)

        # get f(x)
        diff = subVectors(pos, x_v)
        f = dotProduct(diff, diff) - ((c**2) * ((t_v - prev_t)**2))
        # get f'(x)
        f_prime = (2 * dotProduct(subVectors(pos, x_v), getSatelliteVelocity(satellites, index, prev_t))) + (2.0 * (c**2) * (t_v - prev_t))

        cur_t = prev_t - (f / f_prime)

        if abs(cur_t - prev_t) < threshold :
            best_t = cur_t
            break

        prev_t = cur_t

    return getSatelliteLocation(satellites, index, best_t), best_t

# returns a vector that is scaled
def scaleVector(scale, vec) :
    x = vec[0] * scale
    y = vec[1] * scale
    z = vec[2] * scale
    return [x, y, z]

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

# rotates a point along the z-axis
def rotate(u, t_v) :
    time_offset = (2*pi*t_v) / s
    v = [(mp.cos(time_offset) * u[0]) + (-mp.sin(time_offset) * u[1]),
         (mp.sin(time_offset) * u[0]) + (mp.cos(time_offset) * u[1]),
         (u[2])]
    return v

main()