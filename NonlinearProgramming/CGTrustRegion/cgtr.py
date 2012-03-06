#
# Author: Mark Stoehr
#


import numpy as np
from scipy.linalg import norm

def trust_region(x,f,eta,delta_hat,delta,global_tol,max_iter=10000):
    gradient_norm = norm(f.gradient(x))
    num_iterations =0
    while gradient_norm > global_tol and num_iterations <max_iter:
        num_iterations += 1
        if num_iterations % 100 == 0:
            print "In trust region iteration ", num_iterations
            print "Gradient norm is", gradient_norm
            print "rho is", rho
            print "delta is ",delta
            print "m_diff is",m_diff
            print "func_diff is",func_diff
        tol_k = min(.5,np.sqrt(gradient_norm))*gradient_norm
        p,m_diff = steihaug_cg(x,f,tol_k,delta)
        func_diff = f.function(x)-f.function(x+p)
        rho = func_diff/m_diff
        if rho < .25:
            delta = max(.25 * delta,global_tol)
        else:
            if rho  > .75 and (norm(p)-delta) < global_tol:
                delta = min(2*delta,delta_hat)
        if rho > eta:
            x = x + p
            gradient_norm = norm(f.gradient(x))
            
    return x

        
        
    
    

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
    dim = x_k.size
    g_k = f.gradient(x_k)
    B_k = f.hessian(x_k)
    m_k = np.dot(g_k,x_k) \
        +.5 * np.dot(np.dot(B_k,x_k),x_k)
    residual_r = g_k
    direction_d = -g_k
    iterate_z = np.zeros(dim)
    iterate_z_norm = 0.
    num_iterations = 0
    cur_residual_norm = norm(residual_r)
    while True:
        num_iterations += 1
        hess_direction = np.dot(B_k,
                                direction_d)
        d_hess_norm_sq = np.dot(direction_d,
                                hess_direction)
        # check for positive definiteness
        # do line search if not the case
        if d_hess_norm_sq <= 0.:
            p_k,m_p_k = _tau_line_search(iterate_z,
                                    iterate_z_norm,
                                    direction_d,
                                    g_k,
                                    B_k,
                                    delta_k)
            return p_k, m_k -m_p_k

        alpha_j = cur_residual_norm**2/d_hess_norm_sq
        iterate_z_next = iterate_z \
            + alpha_j * direction_d
        
        # check if next iterate goes outside trust region
        iterate_z_next_norm = norm(iterate_z_next)
        if iterate_z_next_norm:
            p_k,m_p_k = _tau_line_search(iterate_z,
                                         iterate_z_norm,
                                         direction_d,
                                         g_k,
                                         B_k,
                                         delta_k)
            return p_k, m_k -m_p_k
        # otherwise update residual
        residual_r_next = residual_r \
            + alpha_j * hess_direction
        next_residual_norm = norm(residual_r_next)
        if next_residual_norm < tol:
            return iterate_z_next

        cur_residual_norm = next_residual_norm
        residual_r = residual_r_next
        direction_d = -residual_r_next \
            +next_residual_norm/cur_residual_norm\
            *direction_d
        iterate_z = iterate_z_next
        iterate_z_norm = iterate_z_next_norm


def _tau_line_search(iterate_z,
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
    if test_value_pos < test_value_neg:
        return test_p_pos,test_value_pos
    else:
        return test_p_neg,test_value_neg
    

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

    


