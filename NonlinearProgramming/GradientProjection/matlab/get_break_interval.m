function [hi,lo,not_found] = get_break_interval(hi,breaks,tol)
  %
  %
  %  In the case where breaks are all greater than tolerance
  %
  %
  %
  %
  lo = hi;
 
  % use exponential search to find hi
  % similar to efficient table expansion technique
  not_found = 1;
  if lo == 0 && breaks(1) > tol
    lo_val = 0;
  else 
    % make sure that lo is at least 1 so that way this is well-defined
    lo = max(1,lo);
    lo_val = breaks(lo);
  end


  num_t = length(breaks);
  % get the boundaries
  hi_lbound = max(lo+1,1);
  % make sure that hi is at least 1 or equal to lo
  hi = hi_lbound;
  hi_ubound = min(lo+20,num_t);

  % search within exponentially growing boundaries (works for infinite streams of values)
  while not_found && (lo <= num_t)
    % check that with respect to the tolerance we distinguish different 
    % entries to make intervals large enough to be interesting
    %display(['Testing with indices ',num2str(hi_lbound),' and ',num2str(hi_ubound)])
   
    % this line is in here to see if it makes it run faster
    hi_ubound = num_t;
    lbound = 0;
    if lo >0
      lbound = breaks(lo);
    end
    ind = find(breaks(hi_lbound:hi_ubound)-lbound > tol + lo_val);
    if  ind 
      hi = ind(1)+hi_lbound-1;
      not_found = 0;
    elseif hi_ubound == num_t
      hi = num_t+1;
      break;
    else
     
      hi_lbound_new = hi_ubound+1;
      hi_ubound = min(num_t,hi_lbound+2*(hi_ubound-hi_lbound));
      hi_lbound = hi_lbound_new;
    end

  end
    
  
