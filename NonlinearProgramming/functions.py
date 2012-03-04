#
#  Author: Mark Stoehr, 2012
#

import numpy as np


# Inner function for the autodiff version
def fenton(x):
    x0 = x[0]
    x1 = x[1]
    y = 12.+x[0]*x[0]\
        +(1.+x[1]*x[1])/(x[0]*x[0])\
        +(x[0]*x[0]*x[1]*x[1]+100.)/((x[0]*x[1])**4)
    y = y/10.
    return y

