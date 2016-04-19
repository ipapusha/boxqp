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
DELTA = 0.1; % added to make KKT system sqd

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
x  = full(x0);
nu = full(nu0);
for iter=1:MAXITER
    % 1. Centering step: minimize t*f0(x) + phi(x) 
    %                    subject to Ax = b
    rd = (qp.d) - (qp.C)*x;
    rp = qp.A*x - qp.b;
    h  = 1./rd;
    H  = [t*(qp.P) + (qp.C)'*spdiags(h.^2,0,nin,nin)*(qp.C), (qp.A)'; (qp.A), sparse(neq, neq)];
    dH = spdiags(DELTA*[ones(nx,1); -ones(neq,1)], 0, nx+neq, nx+neq);
    rhs = -[t*(qp.P*x + qp.q) + (qp.C)'*h; rp];
    
    % solve via LDL factorization
    %sol = H\rhs
    %sol = ldlsparse(H, sd.p, rhs); 
    HpdH = H + dH;
    sol = ldlsparse(HpdH, sd.p, rhs);
    %sol = sol + ldlsparse(HpdH, sd.p, rhs - H*sol); % refinement
end

if k == MAXITER
    fprintf('Failed to converge in MAXITER=%d iterations\n', MAXITER);
end






















end


