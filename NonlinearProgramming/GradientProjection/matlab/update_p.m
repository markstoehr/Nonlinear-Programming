function p=update_p(p,idx)
  %
  % idx - indices corresponding to the ones being set to zero
  % p   - current p iterate
  %
  p(idx) = 0;