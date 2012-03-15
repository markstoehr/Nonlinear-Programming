function tau = get_max_tau (iterate_z,direction_d,u)
  %
  %  We now handle this with the min_tau function
  %   iterate_z + tau * direction_d <= u
  %    
  %   (-iterate_z) + tau * (-direction_d) >= -u
  %
  %
  %
  tau = max( (u - iterate_z)./direction_d);
