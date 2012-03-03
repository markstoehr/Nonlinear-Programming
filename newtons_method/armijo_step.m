function [xx,num_function_evaluations] = armijo_step(fHandle,x,p,rho,beta,g,g_norm,num_function_evaluations)


  f = zeros(100,1);

  f(1) = fHandle(x,0);
  num_function_evaluations = num_function_evaluations+1;
  desc_grad = dot(g, p);
  
  tau = -desc_grad/g_norm^2;
  cur_step = tau * p;

  armijo_constraint = -rho*tau*desc_grad;

  for m=2:100
    xx = x+cur_step;
    f(m) = fHandle(xx,0);
    num_function_evaluations = num_function_evaluations+1;
    armijo_constraint = beta*armijo_constraint;
    % check for the sufficient decrease
    if f(1) - f(m) >= armijo_constraint
      break;
    end
    cur_step = beta *cur_step;
  end