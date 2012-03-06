from numpy import *
import ctypes

lib = ctypes.cdll['./symm_matvec_mult.so']
symm_matvec_mult = lib['symm_matvec_mult']



print "fib(5)", fib(5)

A = array([1,2,3,4])

print "A before doubling", A
print "A's location in memory: %X"%A.ctypes.data
doubleme(A.ctypes.data, len(A))
print "A after doubling", A
