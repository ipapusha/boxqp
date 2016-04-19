% initialize rng
rand('state', 0);
randn('state', 0);

n = 3;
m = 2;
T = 4;

%n = 10;
%m = 3;
%T = 100;

% generate A, B
%A = randn(n,n); 
%B = rand(n,m);
%[U,S,V] = svd(A); 
%A = U*V;
A = reshape(1:n*n,n,n);
B = 10*reshape(1:n*m,n,m);

% generate costs
%Q = eye(n);
%R = eye(m);
Q = reshape(1:n*n,n,n);
Qf = 100*reshape(1:n*n,n,n);
R = reshape((m*m+1):(m*m+m*m),m,m);

% box constraints
xmin = -1*ones(n,1);
xmax = 1*ones(n,1);
umin = -0.5*ones(m,1);
umax = 0.5*ones(m,1);

% initial condition
z = [9+(1:n)]';

% prepare MPC QP
qp = make_mpc(A, B, z, xmin, xmax, umin, umax, Q, Qf, R, T);
sd = analyze_qp(qp);