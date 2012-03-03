function [xx, g_norm,num_function_evaluations,num_linsys_solved] =...
      newton_iter (fHandle, x,rho,beta,num_function_evaluations,num_linsys_solved)
%INPUTS
%  fHandle  - function handle for the one you want to optimize
%             assumed to take two arguments
%  x        - location to start newton step from
%
%OUTPUTS
%  xx       - final value in newton step
%  g_norm   - final gradient at convergence
%

  [f,g,H]=fHandle(x,2);
  num_function_evaluations = num_function_evaluations+1;
  % find the modification tolerance
  H_max_norm = max(sum(abs(H)));
  epsilon = .00001; % much bigger than machine epsilon but ok
  delta = sqrt(epsilon/2)*H_max_norm;

  g_norm = norm(g);

  if g_norm == 0
    xx = x;
  else
   p = ldl_search_direction(H,g,delta);
   p = reshape(p,size(x,1),size(x,2));
  num_linsys_solved = num_linsys_solved+3;
  [xx,num_function_evaluations] = ...
      armijo_step(fHandle,x,p,rho,beta,g,g_norm,num_function_evaluations);
end