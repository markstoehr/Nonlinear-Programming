addpath('../matlab');
addpath('../../../../Optimization/Intlab_V6/');
startintlab
addpath('../../../../Optimization/hw2/');

max_iter = 100000;
tols = [1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,1e-4,1e-3];
problem_sizes = [10,20,40,80,160,320,640,1280,2560,5120,10240];
answers = cell(length(problem_sizes),2);
for k=1:length(problem_sizes)
  n = problem_sizes(k);
  tol = tols(k);
  f=100*ones(n,1);
  l=0.5*ones(n,1);
  u=3.5*ones(n,1);
  x = 2*ones(n,1);
  [~,~,Q] = cute_wrap(x,2);

  [x,num_matvec] = gradient_projection(x,u,l,f,Q,max_iter, tol)
  answers{k,1} = x;
  answers{k,2} = num_matvec;
  save(['up_to_',num2str(k),'.mat'],'answers','x')
end


tols2 = [1e-3,1e-3,1e-3,1e-3,1e-3];
problem_sizes2 = [10240,20480,40960,81920,163840];
answers2 = cell(length(problem_sizes2),2);

for k=1:length(problem_sizes)
  n = problem_sizes2(k);
  tol = tols2(k);
  f=100*ones(n,1);
  l=0.5*ones(n,1);
  u=3.5*ones(n,1);
  x = 2*ones(n,1);
  [~,~,Q] = cute_wrap(x,2);

  [x,num_matvec] = gradient_projection(x,u,l,f,Q,max_iter, tol)
  answers2{k,1} = x;
  answers2{k,2} = num_matvec;
  save(['up_to_',num2str(n),'.mat'],'answers2','x')
end
