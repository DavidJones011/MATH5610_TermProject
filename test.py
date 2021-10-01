import numpy as nm
import sympy as sm
import sys
from decimal import *

def main ():
#    x= sys.stdin.readline()
#    expr = sm.cos(2*x)
#    value = expr.evalf(subs={x:0.0})
#    sys.stdout.write("{}".format(value))

    x = Decimal.from_float(0.72) * Decimal('5')
    sys.stdout.write("{}".format(x))
    pass

main()