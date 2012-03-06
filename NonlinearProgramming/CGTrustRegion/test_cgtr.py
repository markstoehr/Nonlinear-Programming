import cgtr

import os
os.chdir('/home/mark/projects/Nonlinear-Programming')
import NonlinearProgramming as nlpg
from NonlinearProgramming import functions as fun_nlpg
from scipy.linalg import norm

np = fun_nlpg.np
x = fun_nlpg.np.array([3,2])
H_part = np.random.rand(2,2)
quadratic_func = lambda x: np.dot(np.dot(np.dot(H_part,H_part.transpose()),x),x)

H = np.dot(H_part,H_part.transpose())


f = fun_nlpg.Func_Object(quadratic_func,
                         x)

f2 = fun_nlpg.Func_Object(fun_nlpg.fenton,x)

x_k = x
tol = .0001
delta_k = .5
dim = x_k.size
g_k = f.gradient(x_k)
B_k = f.hessian(x_k)
residual_r = g_k
direction_d = -g_k
iterate_z = np.zeros(dim)
iterate_z_norm = 0.
num_iterations = 0
cur_residual_norm = norm(residual_r)


cgtr.steihaug_cg(x_k,f,tol,delta_k)


eta = .1
delta_hat = 10
delta = 1
global_tol = .000001
x_end = cgtr.trust_region(x,f,eta,delta_hat,delta,global_tol,50000)

# homework questions
x,num_iterations,num_matvec, num_function_evaluations = cgtr.trust_region(np.array([3.,2.]),f2,eta,delta_hat,delta,global_tol,50000)

x,num_iterations,num_matvec, num_function_evaluations = cgtr.trust_region(np.array([3.,4.]),f2,eta,delta_hat,delta,global_tol,50000)



x,num_iterations,num_matvec, num_function_evaluations = cgtr.trust_region(np.array([3.,4.]),f2,eta,delta_hat,delta,global_tol,50000,False)

# doing various tests with the cute function
# start with size 10
xs = []
mvs = []
fes = []
for i in [10,100,200,500,1000,5000,10000]:
    x = np.ones(i)
    f = fun_nlpg.Func_Object(fun_nlpg.make_cute_func(x),x)
    x,num_iterations,num_matvec, num_function_evaluations = cgtr.trust_region(x,f,eta,delta_hat,delta,.001,50000,False)
    xs.append(x)
    mvs.append(num_matvec)
    fes.append(num_function_evaluations)

print 

i=10000
x = np.ones(i)
f = fun_nlpg.Func_Object(fun_nlpg.make_cute_func(x),x)
x,num_iterations,num_matvec, num_function_evaluations = cgtr.trust_region(x,f,eta,delta_hat,delta,1,50000,False)


cgtr.steihaug_cg(x,f,.001,.5)

iterate_z = np.random.rand(2)
iterate_z_norm = cgtr.norm(iterate_z)
direction_d = np.random.rand(2)
g_k = f.gradient(x)
B_k = f.hessian(x)
d_hess_norm_sq = np.dot(B_k,direction_d)
delta_k = .5

cgtr.

f2 = 
