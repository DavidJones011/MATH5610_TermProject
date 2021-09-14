# Exercise 12. A "simple" program
# This code is part of an assignmetn for the course MATH5610 taught by Professor Peter Alfeld at the University of Utah
# Author(s): David Jones, Leela, Preston

# This program calcualtes the euclidean distance of a given vector x with n components

import sys
import math

def main() :

    # setup the vector x with n components
    n = int(sys.stdin.readline())

    # calculate the euclidean distance
    distance = 0
    for i in range (0, n) :
        x_i = float(sys.stdin.readline())
        distance += x_i * x_i
    distance = math.sqrt(distance)

    # print result
    sys.stdout.write("{}\n".format(distance))

    pass

main()