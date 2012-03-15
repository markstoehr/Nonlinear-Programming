Q = [ 1 0 
      0 1];
f = [1
     1];

x = [1
     1];
u = [inf
     inf];
l = [-inf
     -inf];
tol = .001;


[x_cauchy,num_matvec_c,A,d,inactive_set,x_inactive] ...
	= get_cauchy_point(x,u,l,f,Q,tol);


t1 = min(x_cauchy == [-1
                      -1]);

G = [ 1 0 
      0 1];
c = [1
     1];

x = [1
     1];
u = [inf
     inf];
l = [-.5
     -.5];
tol = .001;
[x_cauchy,num_matvec_c,A,d,inactive_set,x_inactive] ...
	= get_cauchy_point(x,u,l,c,G,tol);


G = [1 0
     0 1];

c = [0
     0];
x = [10
     10];
u = [11
     11];
l = [2
     2];
