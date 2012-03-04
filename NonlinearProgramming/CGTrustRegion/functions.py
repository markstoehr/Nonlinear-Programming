#
#  Author: Mark Stoehr, 2012
#

import numpy as np


# Inner function for the autodiff version
def _fenton(x):
    x0 = x[0]
    x1 = x[1]
    y = 12.+x[0]*x[0]\
        +(1.+x[1]*x[1])/(x[0]*x[0])\
        +(x[0]*x[0]*x[1]*x[1]+100.)/((x[0]*x[1])**4)
    y = y/10.
    return y




def fenton(x,nderivs):
    """Returns Fenton's function computed at x
    if nderivs is 0 it will return just the function value
    if nderivs is 1 it will return the function value and the
    gradient
    if nderivs is 2 it will return the function value,
    gradient, and the Hessian

    Parameters
    ----------
    x: ndarray
        should be 2 dimensional ndarray
    nderivs: int
        number of derivatives to take

    Output
    ------
    y: double
        function value
    g: ndarray (optional)
        gradient, 2 dimensional ndarray
    h: ndarray (optional)
        hessian, 2 by 2 dimensional ndarray

    """
    y0 = _fenton(x)
    adolc.trace_on(0)
    x = adolc.adouble(x)
    adolc.independent(x)
    y = _fenton(x)
    adolc.dependent(y)
    adolc.trace_off()
    if nderivs == 0:
        return y
    elif nderivs == 1:
        g=adolc.gradient(0,np.ravel(x)).reshape(x.shape)
        return y,g
    elif nderivs == 2:
        g=adolc.gradient(0,np.ravel(x)).reshape(x.shape)
        H=adolc.hessian(0,np.ravel(x)).\
            reshape(x.shape + x.shape)
        return y,g,H
    else:
        print "Error:"
        print "Usage:"
        print " x is 1 dimensonal ndarray of length 2"
        print " nderivs has a value [0,1,2]"

