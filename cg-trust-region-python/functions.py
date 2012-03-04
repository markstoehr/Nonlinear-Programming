#
#  Author: Mark Stoehr, 2012
#


def fenton(x):
    """Returns Fenton's function computed at x
    x should be 2 dimensional ndarray
    """
    x0 = x[0]
    x1 = x[1]
    y = 12.+x1*x1+(1.+x2*x2)/(x1*x1)+(x1*x1*x2*x2+100.)/((x1*x2)^4)
    y = y/10.
    return y


