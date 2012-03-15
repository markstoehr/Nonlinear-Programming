B = eye(4);
g = 4*ones(4,1);

x = ones(4,1);
u = rand(4,1);
u = u/min(u) * 3;
l = rand(4,1);
tol = .001;


[x,num_matvec] = cg_truncated(x,u,l,g,B,tol)

B = eye(4);
g = -4*ones(4,1);

x = ones(4,1);
u = rand(4,1);
u = u/min(u) * 3;
l = rand(4,1);
tol = .001;


[x,num_matvec] = cg_truncated(x,u,l,g,B,tol)
