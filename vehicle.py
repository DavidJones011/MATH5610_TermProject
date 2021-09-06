import sys
import math

from enum import Enum
class GeoCoordinate(Enum) :
    LONGITUDE = 0
    LATITUDE = 1

EARTH_RADIUS = float(6367444.50)
SIDEREAL_DAY_SECS = float(86164.09)

def main() :
    
    for line in sys.stdin :

        if line == '':
            break

        tokens = line.split()

        if(len(tokens) == 0) :
            break

        #store the tokenized data into the correct data types
        t_beg    = float(tokens[0])
        t_end    = float(tokens[1])
        t_steps  = int(tokens[2])
        psi_d    = int(tokens[3])
        psi_m    = int(tokens[4])
        psi_s    = float(tokens[5])
        ns_value = int(tokens[6])
        lambda_d = int(tokens[7])
        lambda_m = int(tokens[8])
        lambda_s = float(tokens[9])
        ew_value = int(tokens[10])
        altitude   = float(tokens[11])

        if(t_steps == 0) :
            t_steps = 1
            delta_time = 1

        if(t_steps > 1) :
            delta_time = (t_end - t_beg) / t_steps

        lat_off = degreeToRad(psi_d, psi_m, psi_s, ns_value, GeoCoordinate.LATITUDE)
        lon_off = degreeToRad(lambda_d, lambda_m, lambda_s, ew_value, GeoCoordinate.LONGITUDE)
        time_rad = 0

        # calculate all positions of the object in each time step
        for step in range(t_steps + 1) :
            time = delta_time * step

            if(t_steps > 1 and step > 0) :
                time_rad = (2.0 * math.pi) * (time / SIDEREAL_DAY_SECS)

            lat_data = radToDegree(lat_off, GeoCoordinate.LATITUDE)
            lon_data = radToDegree(lon_off + time_rad, GeoCoordinate.LONGITUDE)
            sys.stdout.write("{} {} {} {} {} {} {} {} {} {}\n".format(t_beg + t_v, lat_data[0], lat_data[1], lat_data[2], lat_data[3], lon_data[0], lon_data[1], lon_data[2], lon_data[3], altitude))
    pass

def degreeToRad(angle, angle_minute, angle_second, dir_value, geo_coordinate) :
    absolute_angle = angle + (angle_minute / 60.0) + (angle_second / 3600.0)

    # convert the angle to be within (0,360)
    if(geo_coordinate == GeoCoordinate.LONGITUDE) :
        absolute_angle =  ((1 - ((1 + dir_value) * 0.5)) * 360.0) + (dir_value * absolute_angle) 
        
    # convert the angle to be within (0, 180)   
    elif(geo_coordinate == GeoCoordinate.LATITUDE) :
        absolute_angle = 90.0 - (dir_value * absolute_angle)
     
    return absolute_angle * math.pi / 180.0

def radToDegree(angle, geo_coordinate) :
    absolute_angle = angle * 180.0 / math.pi
    dir_value = 1

    # convert the range from (0, 360) to (0, 180)
    # also gets the EW value
    if(geo_coordinate == GeoCoordinate.LONGITUDE) :
        dir_value = int(-1 if (absolute_angle > math.pi) else 1)
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


main()