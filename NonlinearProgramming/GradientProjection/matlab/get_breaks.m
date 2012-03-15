function [breaks,t_bar,t_bar_idx] = get_breaks (x,u,l,g,tol)
  %
  %
  %
  % t_bar records the step size at which point we've hit the boundary for
  % that coordinate
  %
  % breaks is the points that will be in our interval
  % t_bar_idx(k) is the position t_bar for the k-th item
  

  t_bar = max((x-u)./g,(x-l)./g);
  [breaks,t_bar_idx] = sort(t_bar,'ascend');
