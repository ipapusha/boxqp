% initialize rng
rand('state', 0);
randn('state', 0);

%n = 2;
%m = 1;
%T = 4;

n = 8;
m = 4;
T = 20;

% generate A, B
%A = randn(n,n); 
%B = rand(n,m);
%[U,S,V] = svd(A); 
%A = U*V;
A = reshape(1:n*n,n,n)/(n*n);
B = 10*reshape(1:n*m,n,m)/(n*m);

% generate costs
Q = eye(n);
Qf = eye(n);
R = eye(m);
%Q = reshape(1:n*n,n,n);
%Qf = 100*reshape(1:n*n,n,n);
%R = reshape((m*m+1):(m*m+m*m),m,m);

% box constraints
xmin = -1*ones(n,1);
xmax = 1*ones(n,1);
umin = -0.5*ones(m,1);
umax = 0.5*ones(m,1);

% initial condition
z = 0.5*ones(n,1);

% prepare MPC QP
qp = make_mpc(A, B, z, xmin, xmax, umin, umax, Q, Qf, R, T);
ws0 = boxqp_initial(qp);
sd = boxqp_analyze(qp, ws0);
[fval, x, ws, status] = boxqp(qp, sd, ws0);
%[fval, x, ws, status] = boxqp(qp, [], ws0);
fval

%% check timing
REPS = 100; 
tic;
for count=1:REPS
    [fval, x, ws, status] = boxqp(qp, [], ws0, 'quiet');
end
avgtime = toc/REPS;
fprintf('Average solve time (total: %.2f s)\n', avgtime*REPS);
fprintf('  no structure:   %.2f ms (%.1f Hz)\n', ...
        1000*avgtime, 1/avgtime); 

tic;
for count=1:REPS
    [fval, x, ws, status] = boxqp(qp, sd, ws0, 'quiet');
end
avgtime = toc/REPS;
fprintf('  with structure: %.2f ms (%.1f Hz)\n', ...
        1000*avgtime, 1/avgtime);

%% check against CVX 
tic;
cvx_begin;
    variable xx(size(qp.P,1),1);
    minimize( 0.5*quad_form(xx, qp.P) + qp.q'*xx );
        qp.A*xx == qp.b;
        qp.G*xx <= qp.h;
cvx_end;
avgtime = toc/1;
assert(abs(cvx_optval - fval) <= 1e-4);
assert(norm(x-xx) <= 1e-4);
fprintf('        with CVX: %.2f ms (%.1f Hz)\n', ...
        1000*avgtime, 1/avgtime);


%% check actual problem
tic;
cvx_begin;
    variable x(n,T+1);
    variable u(m,T);
    for t=1:T
        x(:,t+1) == A*x(:,t) + B*u(:,t);
        xmin <= x(:,t); x(:,t) <= xmax;
        umin <= u(:,t); u(:,t) <= umax;
    end
    x(:,1) == z;
    xmin <= x(:,T+1); x(:,T+1) <= xmax;
    
    minimize( 0.5*sum(square_pos(norms(sqrtm(Q)*x(:,1:T),2))) + ...
        0.5*sum(square_pos(norms(sqrtm(R)*u,2))) + ...
        0.5*sum(square_pos(norms(sqrtm(Qf)*x(:,T+1),2))) );
cvx_end;
avgtime = toc/1;
fprintf('       naive CVX: %.2f ms (%.1f Hz)\n', ...
        1000*avgtime, 1/avgtime);
