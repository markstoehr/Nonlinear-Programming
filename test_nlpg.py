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


x_cute = np.array([0.158,0.971,0.957,0.485,0.8,0.142,0.422,0.916,0.792,0.959])
cute = fun_nlpg.Func_Object(make_cute_func(x_cute),
                            x_cute)


cute = fun_nlpg.Func_Object(make_cute_func(x_cute),
                            x_cute)
