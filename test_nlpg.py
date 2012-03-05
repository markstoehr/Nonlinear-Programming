import os
os.chdir('/home/mark/projects/Nonlinear-Programming')
import NonlinearProgramming as nlpg
from NonlinearProgramming import functions as fun_nlpg

np = fun_nlpg.np
x = fun_nlpg.np.array([3,2])
f = fun_nlpg.Func_Object(fun_nlpg.fenton,
                         x)

# do an assert equal here
f.function(fun_nlpg.np.array([3,2]))
fun_nlpg.fenton(fun_nlpg.np.array([3,2]))

f.gradient(fun_nlpg.np.array([3,2]))
f.hessian(fun_nlpg.np.array([3,2]))



def make_cute_func(x):
    N = x.shape[0]
    I = np.arange(N-4,dtype=int);
    return lambda z: sum( (-4*z[I]+3.0)**2 ) \
        + sum( ( z[I]**2 + 2*z[I+1]**2\
                   + 3*z[I+2]**2 + 4*z[I+3]**2 \
                   + 5*z[N-1]**2 )**2 );


x_cute = np.ones(100)
cute = fun_nlpg.Func_Object(nlpg.functions.make_cute_func(x_cute),
                            x_cute)

adolc = fun_nlpg.adolc
x_ones_hess = adolc.colpack\
            .sparse_hess_no_repeat(0,np.ravel(np.ones(10)),[0,0])

cute = fun_nlpg.Func_Object(make_cute_func(x_cute),
                            x_cute)

sp_H=cute.sparse_hessian(x_cute)
H=cute.hessian(x_cute)

adolc.colpack\
    .sparse_hess_no_repeat(1,
                           np.ravel(x_cute),
                           [0,0])

adolc.colpack\
    .sparse_hess_repeat(0,
                        np.ravel(x_cute),
                        [0,0])


x_cute2 = np.array([0.158,0.971,0.957,0.485,.23,.5])
cute2 = fun_nlpg.Func_Object(nlpg.functions.make_cute_func(x_cute2),
                            x_cute2)


sp_H=cute2.sparse_hessian(np.ones(6))
H=cute2.hessian(np.ones(6))


x_cute = np.ones(10**4)
cute = fun_nlpg.Func_Object(nlpg.functions.make_cute_func(x_cute),
                            x_cute)

sp_H=cute.sparse_hessian(x_cute)

x_cute = np.ones(10**3*5)
cute = fun_nlpg.Func_Object(nlpg.functions\
                                .make_cute_func(x_cute),
                            x_cute)

sp_H=cute.sparse_hessian(x_cute)




H=cute.hessian(x_cute)
