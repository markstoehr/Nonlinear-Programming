function [x,num_matvec] = tau_line_search(iterate_z,...
				       direction_d,...
				       g,...
				       B,...
				       u,l,...
				       num_matvec)

  % find the optimal value of iterate_z + tau * direction_d to minimize the function and
  % remain within the bounds
  %
  % function to minimize is .5 *(z + tau*d)' * B * (z + tau*d) + (z + tau*d)' * g
  % since we only control tau this is equivalent to minimizing
  %  alpha/2* tau^2  + beta * tau
  %  
  %  so tau = - beta / alpha
  %  
				%  where 
  %    alpha = d'*B*d - positive since B is positive definite and d != 0
  %    beta = d' * B * z + d'g
  %
  %  it needs to be checked whether that is with in the bounds though
  %
  %
  Bd = B * direction_d;
  alpha = direction_d' * Bd;
  zBd = iterate_z' *  Bd;
  dg = direction_d' * g;
  beta = zBd + dg;
  num_matvec = 1 + num_matvec;
  opt_tau = - beta/alpha;
  opt_vec = iterate_z + opt_tau *direction_d;
  if   min( l <= opt_vec) &&  min ( opt_vec <= u) && (alpha >= 0 )
    % last guy is to check that this is a local minimum
    x = opt_vec;
    return
  % otherwise bounds are constraining where we should go
  else 
				% find the smallest tau that gives a feasible point
    bounds = [ (l-iterate_z)./direction_d (u-iterate_z)./direction_d];
    min_tau = max(min(bounds'));
    % find the largest tau that gives a feasible point
    max_tau = min(max(bounds'));

    min_val = .5 * min_tau^2 *alpha + ...
	beta *min_tau;
    max_val = .5 * max_tau^2 *alpha + ...
	beta *max_tau;
    if min_val < max_val
      % display('used min_val')
      x = iterate_z + min_tau*direction_d;
    else
      % display('used max_val')
      x = iterate_z + max_tau*direction_d;
    end
  end