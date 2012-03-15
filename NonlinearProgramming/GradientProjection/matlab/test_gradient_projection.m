Q = eye(4);
f = 4*ones(4,1);

x = ones(4,1);
%u = rand(4,1);
%u = u/min(u) * 3;
u =[    4.3604
    4.3461
    3.0000
    7.9953];

%l = rand(4,1);
l = [  0.3692
    0.1112
    0.7803
    0.3897];


tol = .001;
max_iter = 1000;

[x,num_matvec] = gradient_projection(x,u,l,f,Q,max_iter, tol)

Q = eye(4);
f = -100*ones(4,1);

x = ones(4,1);
%u = rand(4,1);
%u = u/min(u) * 3;
u =[    4.3604
    4.3461
    3.0000
    7.9953];

%l = rand(4,1);
l = [  0.3692
    0.1112
    0.7803
    0.3897];


tol = .000001;
max_iter = 1000;

[x,num_matvec] = gradient_projection(x,u,l,f,Q,max_iter, tol)
