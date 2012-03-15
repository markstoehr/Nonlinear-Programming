function tau = get_min_tau (iterate_z,direction_d,l)

  %  We want
  %    iterate_z + tau * direction_d >= l
  %    
  %    tau* direction_d >= l- iterate_z
  %    
  %    tau = max( (l- iterate_z) ./ direction_d)
  %
  %
  tau = min( (l - iterate_z)./direction_d);