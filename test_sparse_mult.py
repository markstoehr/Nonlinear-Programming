import os
os.chdir('/home/mark/projects/Nonlinear-Programming')
import NonlinearProgramming as nlpg
from NonlinearProgramming import functions as fun_nlpg
from numpy import triu

np = fun_nlpg.np


def to_sp_hess(hess):
    hess_triu = triu(H)
    nnzs = hess_triu.nonzero()
    num_idx = nnzs[0].shape
    return [num_idx[0],nnzs[0],nnzs[1],H[nnzs]]

# test that the sparse function works
# did a test using other code

x_cute = np.ones(7)
cute = fun_nlpg.Func_Object(nlpg.functions.make_cute_func(x_cute),
                            x_cute)


sp_H=cute.sparse_hessian(x_cute)
H=cute.hessian(x_cute)

##
#  Need an array testing scheme right here
#
#
