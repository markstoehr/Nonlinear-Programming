function [x, g_norms,num_function_evaluations,num_linsys_solved] = newton_method (fHandle, x_0, tol, max_iter,rho,beta)
%INPUTS
%  fHandle  - function handle for the one you want to optimize
%             assumed to take two arguments
%  x_0      - start point
%  tol      - max value of gradient of function before we terminate
%  max_iter - max number of iterations
%  rho/beta  - parameters for the armijo line search rule
%OUTPUTS
%  xx       - final value at convergence
%  g_norm   - final gradient at convergence
%

assert (tol > 0 ...
	&& max_iter > 1);

g_norms = zeros(max_iter,1);

% current iterate
x = x_0;

do_message = 1;

num_function_evaluations = 0;
num_linsys_solved = 0;
for  k=1:max_iter
  [x,g_norm,num_function_evaluations,num_linsys_solved] = ...
      newton_iter (fHandle, x,rho,beta,num_function_evaluations,num_linsys_solved);
  g_norms(k) = g_norm;
  if g_norm < tol
    break;
  end
end

g_norms = g_norms(1:k);
display(['Optimization ended after ',num2str(k),' iterations with norm of gradient = ',num2str(g_norm)]);
display(['Number of function evaluations: ',num2str(num_function_evaluations)]);
display(['Number of linear systems solved: ',num2str(num_linsys_solved)]);
