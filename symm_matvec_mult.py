from numpy import *
import ctypes
import NonlinearProgramming as nlpg
from NonlinearProgramming import functions as fun_nlpg

np = fun_nlpg.np

x_cute = np.ones(7)
cute = fun_nlpg\
    .Func_Object(nlpg.functions\
                     .make_cute_func(x_cute),
                 x_cute)


sp_H=cute.sparse_hessian(x_cute)
H=cute.hessian(x_cute)



lib = ctypes.cdll['./symm_matvec_mult.so']
symm_matvec_mult = lib['symm_matvec_mult']

Hridx = sp_H[1]
Hcidx = sp_H[2]
Hval = sp_H[3]
nIdx = sp_H[0]

inVec = np.random.rand(7)
outVec = np.zeros(7)


print "Hridx's location in memory: %X"%Hridx.ctypes.data
print "Hcidx's location in memory: %X"%Hcidx.ctypes.data
print "Hval's location in memory: %X"%Hval.ctypes.data
print "inVec's location in memory: %X"%inVec.ctypes.data
print "outVec's location in memory: %X"%outVec.ctypes.data

print "outVec before filling", outVec

symm_matvec_mult(Hridx.ctypes.data,
                 Hcidx.ctypes.data,
                 Hval.ctypes.data,
                 inVec.ctypes.data,
                 outVec.ctypes.data,
                 nIdx)

print "outVec after filling", outVec


print "fib(5)", fib(5)

A = array([1,2,3,4])

print "A before doubling", A
print "A's location in memory: %X"%A.ctypes.data
doubleme(A.ctypes.data, len(A))
print "A after doubling", A
