#
# Author: Mark Stoehr 2012
#
#

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
        y = _fenton(x)
        adolc.dependent(y)
        adolc.trace_off()
        self.tape_number = CUR_ADOLC_TAPE_NUMBER
        CUR_ADOLC_TAPE_NUMBER += 1

    def function(x):
        return adolc.function(self.tape_number,
                              np.ravel(x))[0]
    
    def gradient(x):
        return adolc.gradient(self.tape_number,
                              np.ravel(x)).reshape(x.shape)
    
    def hessian(x):
        return adolc.hessian(self.tape_number,
                             np.ravel(x)).reshape(x.shape)


