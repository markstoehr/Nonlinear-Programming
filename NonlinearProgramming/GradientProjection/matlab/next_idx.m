function right_idx = next_idx(breaks,left_idx,tol)

  right_idx = left_idx+1;
  while abs(breaks(right_idx) - breaks(left_idx)) < tol
    right_idx = right_idx +1;
  end