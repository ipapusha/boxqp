function [fval,x] = box_qp(qp, sd, x0, nu0)
% extract dimensions
nx  = size(qp.P, 1); % number of variables
neq = size(qp.A, 1); % number of equality constraints
nin = size(qp.C, 1); % number of inequality constraints


% algorithm parameters
MU  = 10;
TOL = 1e-4;
ALPHA = 0.01;
BETA  = 0.5;
MAXITER = 100;

% Barrier method
% given: strictly feasible x, 
%        t > 0, mu > 1, tolerance epsilon > 0
% repeat:
% 1. Centering step: 
%    compute x^*(t) by minimizing t*f0(x) + phi(x) subject to Ax = b, starting at x
% 2. Update: x := x^*(t)
% 3. Stopping criterion: quit if m/t <= eps
% 4. Increase t: t = mu*t

% Infeasible start Newton's method (Algorithm 10.2)
% given: starting point x, nu, 
%        tolerance epsilon > 0, 
%        0 < alpha < 1/2, 
%        0 < beta < 1
% repeat:
% 1. Compute primal and dual Newton steps 
% 2. Backtracking line search 
% 3. Update

t  = 1.1;
x  = sparse(x0);
nu = sparse(nu0);
for iter=1:MAXITER
    % 1. Centering step: minimize t*f0(x) + phi(x) 
    %                    subject to Ax = b
    rd = (qp.d) - (qp.C)*x;
    rp = qp.A*x - qp.b;
    h = 1./rd;
    H = [t*(qp.P) + (qp.C)'*spdiags(h.^2,0,nin,nin)*(qp.C), (qp.A)'; (qp.A), sparse(neq, neq)];
    rhs = -full([t*(qp.P*x + qp.q) + (qp.C)'*h; rp]);
    sol = ldlsparse(H, sd.p, rhs); % sol = H\rhs
end

if k == MAXITER
    fprintf('Failed to converge in MAXITER=%d iterations\n', MAXITER);
end




























end

