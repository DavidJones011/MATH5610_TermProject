import numpy as nm
import sympy as sm
import sys
from decimal import *
import mpmath as mp

def main ():
#    x= sys.stdin.readline()
#    expr = sm.cos(2*x)
#    value = expr.evalf(subs={x:0.0})
#    sys.stdout.write("{}".format(value))

    mp.mp.prec = 20

    x = mp.mpf(55.0 / 3600.0)
    y = mp.mpf(55 / 3600.0)
    sys.stdout.write("{}, {}".format(x, y))
    pass

main()