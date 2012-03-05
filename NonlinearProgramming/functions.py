#
#  Author: Mark Stoehr, 2012
#

import adolc
import numpy as np
from scipy.sparse import coo_matrix

# constant for which tape number to use
CUR_ADOLC_TAPE_NUMBER = 0

class Func_Object:
    def __init__(self,func,x):
        """ Initializes the tape for a function which allows
        the user to then quickly evaluate the function,
        gradient, and hessian at arbitrary points

        func is assumed to be a function from a vector array
        to a single number.

        Also updates the global CUR_ADOLC_TAPE_NUMBER which
        keeps track of the tapes that have been used by
        ADOLC

        Parameters
        ----------
        func: function
            function to be evaluated
        x: ndarray
            ndarray of correct shape and type for the 
            function

        """
        # access the tape number counter
        global CUR_ADOLC_TAPE_NUMBER
        # initialize the ADOLC tape for the function
        adolc.trace_on(CUR_ADOLC_TAPE_NUMBER)
        x = adolc.adouble(x)
        adolc.independent(x)
        y = func(x)
        adolc.dependent(y)
        adolc.trace_off()
        self.tape_number = CUR_ADOLC_TAPE_NUMBER
        CUR_ADOLC_TAPE_NUMBER += 1

    def function(self,x):
        return adolc.function(self.tape_number,
                              np.ravel(x))[0]
    
    def gradient(self,x):
        return adolc.gradient(self.tape_number,
                              np.ravel(x)).reshape(x.shape)
    
    def hessian(self,x):
        return adolc.hessian(self.tape_number,
                             np.ravel(x))\
                             .reshape(x.shape\
                                          + x.shape)
    
    def sparse_hessian(self,x):
        sp_hess = adolc.colpack\
            .sparse_hess_no_repeat(self.tape_number,
                                   np.ravel(x),
                                   [0,0])
        # sparse matrices are stored in csc format
        # in order to do fast arithmetic
        return sp_hess
    #coo_matrix(sp_hess[3],
     #                     (sp_hess[1],
     #                      sp_hess[2])).tocsc()
            


class Func_Object_Sparse:
    def __init__(self,func,x):
        """ Initializes the tape for a function which allows
        the user to then quickly evaluate the function,
        gradient, and hessian at arbitrary points with
        a sparsity pattern

        func is assumed to be a function from a vector array
        to a single number.

        Also updates the global CUR_ADOLC_TAPE_NUMBER which
        keeps track of the tapes that have been used by
        ADOLC

        Parameters
        ----------
        func: function
            function to be evaluated
        x: ndarray
            ndarray of correct shape and type for the 
            function

        """
        # access the tape number counter
        global CUR_ADOLC_TAPE_NUMBER
        # initialize the ADOLC tape for the function
        adolc.trace_on(CUR_ADOLC_TAPE_NUMBER)
        x = adolc.adouble(x)
        adolc.independent(x)
        y = func(x)
        adolc.dependent(y)
        adolc.trace_off()
        self.tape_number = CUR_ADOLC_TAPE_NUMBER
        CUR_ADOLC_TAPE_NUMBER += 1

    def function(self,x):
        return adolc.function(self.tape_number,
                              np.ravel(x))[0]
    
    def gradient(self,x):
        return adolc.gradient(self.tape_number,
                              np.ravel(x)).reshape(x.shape)
    
    def hessian(self,x):
        return adolc.colpack.sparhessian(self.tape_number,
                             np.ravel(x)).reshape(x.shape\
                                                      + x.shape)


# Inner function for the autodiff version
def fenton(x):
    x0 = x[0]
    x1 = x[1]
    y = 12.+x[0]*x[0]\
        +(1.+x[1]*x[1])/(x[0]*x[0])\
        +(x[0]*x[0]*x[1]*x[1]+100.)/((x[0]*x[1])**4)
    y = y/10.
    return y

def make_cute_func(x):
    N = x.shape[0]
    I = np.arange(N-4,dtype=int);
    return lambda z: sum( (-4*z[I]+3.0)**2 ) \
        + sum( ( z[I]**2 + 2*z[I+1]**2\
                   + 3*z[I+2]**2 + 4*z[I+3]**2 \
                   + 5*z[N-1]**2 )**2 );
