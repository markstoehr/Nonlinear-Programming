function [x,num_matvec] = cg_truncated(x,u,l,g,B,tol)
  % INPUT:
  %
  %  We aim to minimize .5 *x' * B * x + x' * g    s.t. l <= x <= u
  %
  %
  %  x       - vector
  %  g       - vector
  %  B       - matrix
  %  tol     - tolerance for stopping
  %  u,l     - upper and lower bounds, respectively
  %
  % OUTPUT:
  %  x          - vector
  %  num_matvec - number of matrix-vector multiplications
  %
  %

  dim = length(x);
  if dim == 0
    num_matvec = 0;
    return
  end

  residual_r = B*x + g;
  norm_gradient = norm(residual_r);
  direction_d = -residual_r;
  iterate_z = x;
  iterate_z_next = x;

  num_iterations = 0;
  num_matvec = 0;
  cur_residual_norm_sq = residual_r' * residual_r;

  while norm_gradient > tol
    % Hessian calculations
    hess_d = B*direction_d;
    d_hess_norm_sq = direction_d' * hess_d;
    num_matvec = num_matvec+1;

    % exit if the norm of the gradient
    % as measured by the hessian is small
    if d_hess_norm_sq < tol
      [x,num_matvec] = tau_line_search(iterate_z, ...
				       direction_d,...
				       g,...
				       B,...
				       u,l,...
				       num_matvec);
      return;
    end

    alpha= cur_residual_norm_sq/d_hess_norm_sq;
    iterate_z_next = iterate_z + alpha *direction_d;
    

    % we see if we have gone beyond the boundary
    % with the step induced by alpha
    if ~min(min (iterate_z_next <= u), min(iterate_z_next>=l))
      [x,num_matvec] = tau_line_search(iterate_z, ...
				       direction_d,...
				       g,...
				       B,...
				       u,l,...
				       num_matvec);
      return;
      
    end
    
    residual_r_next = residual_r + alpha *hess_direction;
    next_residual_norm_sq = residual_r_next' * residual_r_next;
    
    if next_residual_norm_sq < tol^2
      x = iterate_z_next;
      return;
    end
    

    residual_r = residual_r_next;
    direction_d = -residual_r+...
	(next_residual_norm_sq/cur_residual_norm_sq)*direction_d;
    cur_residual_norm_sq = next_residual_norm_sq;
    iterate_z = iterate_z_next;
    norm_gradient = norm(B * iterate_z + g);
    num_matvec = num_matvec+1;
  end 
  x = iterate_z_next;