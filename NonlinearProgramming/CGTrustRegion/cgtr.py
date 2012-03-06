#
# Author: Mark Stoehr
#


import numpy as np
from scipy.linalg import norm

def trust_region(x,f,eta,delta_hat,delta,global_tol,max_iter=10000,useDense=True):
    gradient_norm = norm(f.gradient(x))
    num_iterations =0
    num_matvec = 0
    num_function_evaluations = 0
    f_x = f.function(x)
    while gradient_norm > global_tol and num_iterations <max_iter:
        num_iterations += 1
        if num_iterations % 100 == 0:
            print "In trust region iteration ", num_iterations
            print "Gradient norm is", gradient_norm
            print "rho is", rho
            print "delta is ",delta
            print "m_diff is",m_diff
            print "func_diff is",func_diff
            print "p is", p
        tol_k = min(.5,np.sqrt(gradient_norm))*gradient_norm
        print "doing cg iteration"
        if useDense:
            p,m_diff,add_matvec_mult = steihaug_cg(x,f,tol_k,delta)
        else:
            p,m_diff,add_matvec_mult = sp_steihaug_cg(x,f,tol_k,delta)
        num_matvec += add_matvec_mult
        f_x_new = f.function(x+p)
        num_function_evaluations +=1
        func_diff = f_x-f_x_new
        rho = func_diff/m_diff
        print "checking rho"
        if rho < .25:
            print "rho<.25"
            delta = max(.25 * delta,global_tol)
        else:
            print "rho>.25"
            if rho  > .75 and (norm(p)-delta) < global_tol:
                delta = min(2*delta,delta_hat)
        print "checking eta"
        if rho > eta:
            print "rho>eta"
            x = x + p
            f_x = f_x_new
            gradient_norm = norm(f.gradient(x))
            print "updated gradient norm"
    return x,num_iterations,num_matvec, num_function_evaluations

# sparse steihaug iteration        
def sp_steihaug_cg(x_k,f,tol,delta_k):
    """Use Steihaug's conjugated gradients algorithm to
    minimize the trust region objective

    Parameters
    ----------
    x_k: array
        initializing point
    f: func_object
        Hessian (or approximation) at the point
    tol: float64
        tolerance 
    """
    print "tol for cg is", tol
    dim = x_k.size
    g_k = f.gradient(x_k)
    sp_B_k = f.sparse_hessian(x_k)
    num_matvec=0
    residual_r = g_k
    direction_d = -g_k
    iterate_z = np.zeros(dim)
    iterate_z_norm = 0.
    num_iterations = 0
    cur_residual_norm = norm(residual_r)
    for i in xrange(20):
        if num_iterations % 50 == 0:
            print "on iteration %d of steihaug cg" % num_iterations
        num_iterations += 1
        hess_direction = f.sparse_symm_matvec_mult(sp_B_k,
                                direction_d)
        d_hess_norm_sq = np.dot(direction_d,
                                hess_direction)
        num_matvec+=1
        # check for positive definiteness
        # do line search if not the case
        print "checking d_hess_norm"
        if d_hess_norm_sq <= 0.:
            return  _sp_tau_line_search(iterate_z,
                                        iterate_z_norm,
                                        direction_d,
                                        g_k,
                                        sp_B_k,
                                        delta_k,num_matvec,
                                        f.sparse_symm_matvec_mult)
            

        alpha_j = cur_residual_norm**2/d_hess_norm_sq
        iterate_z_next = iterate_z \
            + alpha_j * direction_d
        
        # check if next iterate goes outside trust region
        iterate_z_next_norm = norm(iterate_z_next)
        if iterate_z_next_norm >= delta_k:
            return  _sp_tau_line_search(iterate_z,
                                        iterate_z_norm,
                                        direction_d,
                                        g_k,
                                        sp_B_k,
                                        delta_k,
                                        num_matvec,
                                        f.sparse_symm_matvec_mult)

        # otherwise update residual
        residual_r_next = residual_r \
            + alpha_j * hess_direction
        next_residual_norm = norm(residual_r_next)
        print "residual norm is ", next_residual_norm
        if next_residual_norm < tol:
            # our returned p is now iterate_z_next
            # we also predict the reduction
            neg_m_diff =  np.dot(g_k,iterate_z_next) + .5\
                *np.dot(\
                f.sparse_symm_matvec_mult(sp_B_k,
                                          iterate_z_next),
                iterate_z_next)
            num_matvec+=1    
            return iterate_z_next, -neg_m_diff, num_matvec

        cur_residual_norm = next_residual_norm
        residual_r = residual_r_next
        direction_d = -residual_r_next \
            +next_residual_norm/cur_residual_norm\
            *direction_d
        iterate_z = iterate_z_next
        iterate_z_norm = iterate_z_next_norm
    neg_m_diff = np.dot(g_k,iterate_z_next) + .5\
        *np.dot(\
        f.sparse_symm_matvec_mult(sp_B_k,iterate_z_next),
        iterate_z_next)
    num_matvec +=1
    return iterate_z_next, -neg_m_diff,num_matvec


def _sp_tau_line_search(iterate_z,
                     iterate_z_norm,
                     direction_d,
                     g_k,
                     sp_B_k,
                     delta_k,num_matvec,sparse_mult_func):
    """Find a tau that minimizes
    <(g_k),(z+tau*d)> +1/2*||(z+tau*d)||_{sp_B_k}

    """
    dz_prod = np.dot(iterate_z,direction_d)
    d_norm_sq = np.dot(direction_d,direction_d)
    z_norm_sq = iterate_z_norm**2
    discriminant = np.sqrt(dz_prod**2 \
                               -d_norm_sq*(z_norm_sq\
                                               -delta_k**2))\
                                               /d_norm_sq
    center = -dz_prod/d_norm_sq
    tau_pos = center+discriminant
    tau_neg = center-discriminant
    # solution is on either side of the trust region along
    # the line describe by z + tau * d
    test_p_pos = iterate_z+tau_pos*direction_d
    test_p_neg = iterate_z+tau_neg*direction_d
    # find which end of the trust region gives the best
    # result
    test_value_pos = np.dot(g_k,test_p_pos)\
        +.5*np.dot(\
        sparse_mult_func(sp_B_k,
                         test_p_pos),
        test_p_pos)
    test_value_neg = np.dot(g_k,test_p_neg)\
        +.5*np.dot(\
        sparse_mult_func(sp_B_k,
                         test_p_neg),
        test_p_neg)
    num_matvec +=2
    if test_value_pos < test_value_neg:
        return test_p_pos,-test_value_pos,num_matvec
    else:
        return test_p_neg,-test_value_neg,num_matvec
    
        
    
    

def steihaug_cg(x_k,f,tol,delta_k):
    """Use Steihaug's conjugated gradients algorithm to
    minimize the trust region objective

    Parameters
    ----------
    x_k: array
        initializing point
    f: func_object
        Hessian (or approximation) at the point
    tol: float64
        tolerance 
    """
    print "tol for cg is", tol
    dim = x_k.size
    g_k = f.gradient(x_k)
    B_k = f.hessian(x_k)
    num_matvec=0
    residual_r = g_k
    direction_d = -g_k
    iterate_z = np.zeros(dim)
    iterate_z_norm = 0.
    num_iterations = 0
    cur_residual_norm = norm(residual_r)
    for i in xrange(20):
        if num_iterations % 50 == 0:
            print "on iteration %d of steihaug cg" % num_iterations
        num_iterations += 1
        hess_direction = np.dot(B_k,
                                direction_d)
        d_hess_norm_sq = np.dot(direction_d,
                                hess_direction)
        num_matvec+=1
        # check for positive definiteness
        # do line search if not the case
        print "checking d_hess_norm"
        if d_hess_norm_sq <= 0.:
            return  _tau_line_search(iterate_z,
                                    iterate_z_norm,
                                    direction_d,
                                    g_k,
                                    B_k,
                                    delta_k,num_matvec)
            

        alpha_j = cur_residual_norm**2/d_hess_norm_sq
        iterate_z_next = iterate_z \
            + alpha_j * direction_d
        
        # check if next iterate goes outside trust region
        iterate_z_next_norm = norm(iterate_z_next)
        if iterate_z_next_norm >= delta_k:
            return  _tau_line_search(iterate_z,
                                     iterate_z_norm,
                                     direction_d,
                                     g_k,
                                     B_k,
                                     delta_k,num_matvec)

        # otherwise update residual
        residual_r_next = residual_r \
            + alpha_j * hess_direction
        next_residual_norm = norm(residual_r_next)
        print "residual norm is ", next_residual_norm
        if next_residual_norm < tol:
            # our returned p is now iterate_z_next
            # we also predict the reduction
            neg_m_diff =  np.dot(g_k,iterate_z_next) + .5\
                *np.dot(np.dot(B_k,iterate_z_next),
                        iterate_z_next)
            num_matvec+=1    
            return iterate_z_next, -neg_m_diff, num_matvec

        cur_residual_norm = next_residual_norm
        residual_r = residual_r_next
        direction_d = -residual_r_next \
            +next_residual_norm/cur_residual_norm\
            *direction_d
        iterate_z = iterate_z_next
        iterate_z_norm = iterate_z_next_norm
    neg_m_diff = np.dot(g_k,iterate_z_next) + .5\
                *np.dot(np.dot(B_k,iterate_z_next),
                        iterate_z_next)
    num_matvec +=1
    return iterate_z_next, -neg_m_diff,num_matvec


def _tau_line_search(iterate_z,
                     iterate_z_norm,
                     direction_d,
                     g_k,
                     B_k,
                     delta_k,num_matvec):
    """Find a tau that minimizes
    <(g_k),(z+tau*d)> +1/2*||(z+tau*d)||_{B_k}

    """
    dz_prod = np.dot(iterate_z,direction_d)
    d_norm_sq = np.dot(direction_d,direction_d)
    z_norm_sq = iterate_z_norm**2
    discriminant = np.sqrt(dz_prod**2 \
                               -d_norm_sq*(z_norm_sq\
                                               -delta_k**2))\
                                               /d_norm_sq
    center = -dz_prod/d_norm_sq
    tau_pos = center+discriminant
    tau_neg = center-discriminant
    test_p_pos = iterate_z+tau_pos*direction_d
    test_p_neg = iterate_z+tau_neg*direction_d
    test_value_pos = np.dot(g_k,test_p_pos)\
        +.5*np.dot(np.dot(B_k,test_p_pos),test_p_pos)
    test_value_neg = np.dot(g_k,test_p_neg)\
        +.5*np.dot(np.dot(B_k,test_p_neg),test_p_neg)
    num_matvec +=2
    if test_value_pos < test_value_neg:
        return test_p_pos,-test_value_pos,num_matvec
    else:
        return test_p_neg,-test_value_neg,num_matvec
    

def _tau_line_search_diagnose(iterate_z,
                     iterate_z_norm,
                     direction_d,
                     g_k,
                     B_k,
                     delta_k):
    """Find a tau that minimizes
    <(g_k),(z+tau*d)> +1/2*||(z+tau*d)||_{B_k}

    """
    dz_prod = np.dot(iterate_z,direction_d)
    d_norm_sq = np.dot(direction_d,direction_d)
    z_norm_sq = iterate_z_norm**2
    discriminant = np.sqrt(dz_prod**2 \
                               -d_norm_sq*(z_norm_sq\
                                               -delta_k**2))\
                                               /d_norm_sq
    center = -dz_prod/d_norm_sq
    tau_pos = center+discriminant
    tau_neg = center-discriminant
    test_p_pos = iterate_z+tau_pos*direction_d
    test_p_neg = iterate_z+tau_neg*direction_d
    test_value_pos = np.dot(g_k,test_p_pos)\
        +.5*np.dot(np.dot(B_k,test_p_pos),test_p_pos)
    test_value_neg = np.dot(g_k,test_p_neg)\
        +.5*np.dot(np.dot(B_k,test_p_neg),test_p_neg)
    return test_p_pos, test_p_neg

    


