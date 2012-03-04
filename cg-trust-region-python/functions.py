#
#  Author: Mark Stoehr, 2012
#



# Inner function for the autodiff version
def _fenton(x):
    x0 = x[0]
    x1 = x[1]
    y = 12.+x1*x1+(1.+x2*x2)/(x1*x1)+(x1*x1*x2*x2+100.)/((x1*x2)^4)
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
    

