function [x,num_matvec] = gradient_projection(x,u,l,f,Q,max_iter, tol)


  num_matvec = 1;
  [satisfied,max_violation] = check_KKT(Q*x + f,u,l,x, tol);
  for k=1:max_iter
    if satisfied
      display(['Optimization completed after ',num2str(k),' iterations '])
      return
    end
    if mod(k-1,500)==0
      display(['On iteration ',num2str(k),' with max_violation of ',num2str(max_violation)]);
      
    end
    [x_cauchy,num_matvec_c,A,d,inactive_set,x_inactive] ...
	= get_cauchy_point(x,u,l,f,Q,tol);
    num_matvec = num_matvec_c+num_matvec;
    [x_inactive,num_matvec_t] ...
	= cg_truncated(x_inactive,u(inactive_set),...
		       l(inactive_set),d,A,tol);
    x = x_cauchy;
    x(inactive_set) = x_inactive;
    num_matvec = num_matvec_t + num_matvec;

    [satisfied,max_violation] = check_KKT(Q*x + f,u,l,x, tol);
    num_matvec = num_matvec +1;
  end
    

  