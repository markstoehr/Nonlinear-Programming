function f_prime = get_f_prime(c,G,p,x)
  f_prime = c'*p + x' * G * p;