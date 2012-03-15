function [x,num_matvec,A,d,inactive_set,x_inactive] = get_cauchy_point(x,u,l,c,G,tol)

  num_matvec = 0;
  %gradient is found using the two parameters for the quadratic objective
  %display(['G is ',num2str(size(G))])
  g = G*x+c;

  num_matvec = num_matvec +1;
  % make sure that we can make a nontrivial step
  if norm(g) < tol
    % gradient is arbitrarily small
    A = zeros(0,0);
    d = zeros(0,1);
    inactive_set = zeros(0,1);
    x_inactive = zeros(0,1);
    return;
  end
  
  % breaks are the points where if we take a step in the direction with g we will hit
  % boundary
  [breaks,t_bar,t_bar_idx] = get_breaks(x,u,l,g);
  
  % we have two cases

  % find where the interval indices are
  % hi is the smallest index which has greater than tolerance value of lo
  % function implicitly begins with a shadow previous iterate where hi = 1
  [hi,lo] = get_break_interval(0,breaks,min(tol,1/ max(abs(g))));
  if hi == 0
    display(['hi is oddly equal to 0'])
  end

  % we initialize p to be -g, and we set the components corresponding
  %display(['bounds are ',num2str(lo),' and ',num2str(hi-1)])
  if lo == 0
    p=update_p(-g,1:0);
  else
    p=update_p(-g,t_bar_idx(lo:hi-1));
  end

  f_prime = get_f_prime(c,G,p,x);
  f_dprime = get_f_dprime(p,G);
  num_matvec = num_matvec +2;
  gamma = - f_prime / f_dprime;
  % search for the cauchy
  cauchy_not_found =1;
  num_iterations = 0;
  % get x(0) is equal just to x

  while 1
    if mod(hi,20)==0
      display(['Number of intervals processed has been ',num2str(hi)])
      display(['f_prime is ',num2str(f_prime)])
    end
    num_iterations = num_iterations +1;
    if lo == 0
      t_j = 0;
    else
      t_j = breaks(lo);
    end
    t_jp1 = breaks(min(end,hi));
    
    % assume that x = x(t_j)
    if f_prime > -tol
      x = x+ t_j*p;
      break;
    elseif ((0 <= gamma) && (gamma < (t_jp1 - t_j)))
      x = x+gamma*p;
      break;
    elseif            ((t_jp1-t_j) < tol)
      % already updated x at the end of the previous while loop iteration
      break;
    end

				% get the two bounds on the interval
    
    

    
    % update x to x(t_jp1) for the next iteration
    x = x + (t_jp1-t_j)*p;
    [hi,lo] = get_break_interval(hi,breaks,min(tol,1/ max(abs(g))));
    p=update_p(p,t_bar_idx(lo:hi-1));
    f_prime = get_f_prime(c,G,p,x);
    f_dprime = get_f_dprime(p,G);
    num_matvec = num_matvec +2;
    gamma = - f_prime / f_dprime;

  end

  
  inactive_set = find(p);

  a = ones(size(p));
  a(inactive_set) = 0;
  active_set = find(a);
  
  A = G(inactive_set,inactive_set);
  C = G(inactive_set,active_set);
  B = G(active_set,active_set);
  
  % represents new direction
  y = x(active_set);
  d = C*y + c(inactive_set);
  x_inactive = x(inactive_set);
  