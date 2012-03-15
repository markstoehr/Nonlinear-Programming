function [satisfied, max_violation] = check_KKT(grad_func, u,l,x,tol)
  at_ubound = abs(x-u) < tol;
  at_lbound = abs(x-l) < tol;
  between_bounds = min(1-at_ubound, 1-at_lbound);

  idx_at_ubound= find( at_ubound);
  idx_at_lbound = find( at_lbound);
  idx_between = find(between_bounds);

  ubound_test = min(grad_func(idx_at_ubound) < tol);
  if length(ubound_test) == 0
    ubound_test = 1;
  end
  lbound_test = min(grad_func(idx_at_lbound) > -tol);
  if length(lbound_test) == 0
    lbound_test = 1;
  end
  between_test = min(abs(grad_func(idx_between)) < tol);
  if length(between_test) == 0
    between_test = 1;
  end

  satisfied = min(min(between_test,ubound_test),lbound_test);
  max_violation = -inf;
  if ~satisfied
    if ~ubound_test
      max_violation = max(max_violation,max(grad_func(idx_at_ubound)));
    end
    if ~lbound_test
      max_violation = max(max_violation,max(-grad_func(idx_at_ubound)));
    end
    if ~between_test
      max_violation = max(max_violation,max(abs(grad_func(idx_between))));
    end
  end