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
DELTA = 0.1;   % added to make KKT system sqd

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
while neq/t > TOL
    % Centering: minimize t*f0(x) + phi(x) 
    %            subject to Ax = b
    fprintf('Centering step\n');
    
    for iter=1:MAXITER
        % 1. compute primal and dual Newton steps
        z  = qp.d - qp.C*x; h = 1./z;
        rp = qp.A*x - qp.b;
        rd = t*(qp.P*x + qp.q) + qp.C'*h + qp.A'*nu;
        pdgap = norm([rp;rd]);
        
        if pdgap <= TOL
            break;
        end
        
        H   = [t*(qp.P) + qp.C'*spdiags(h.^2,0,nin,nin)*qp.C, qp.A'; qp.A, sparse(neq, neq)];
        dH  = [DELTA*speye(nx), sparse(nx,neq); sparse(neq,nx), -DELTA*speye(neq)];
        rhs = -[rd; rp];

        sol = H\rhs;
        %HpdH = H + dH;
        %sol  = ldlsparse(HpdH, sd.p, rhs);
        %sol = sol + ldlsparse(HpdH, sd.p, rhs - H*sol); % refinement step
        dx = sol(1:nx); dnu = sol(nx+1:nx+neq);

        % 2. Backtracking line search
        tstep = 1;
        while 1
            % compute r(x+t*dx, nu+t*dnu)
            z1  = qp.d - qp.C*(x+tstep*dx); h1 = 1./z1;
            rp1 = qp.A*(x+tstep*dx) - qp.b;
            rd1 = t*(qp.P*(x+tstep*dx) + qp.q) + qp.C'*h1 + qp.A'*(nu+tstep*dnu);

            % check if step size leads to suitable gap decrement
            if norm([rp1; rd1], 2) <= (1-ALPHA*tstep)*pdgap
                break;
            end

            % decrease step
            tstep = BETA*tstep;
        end
        
        fprintf('iter: %3d rp: %6g, rd: %6g, gap: %6g, step: %g\n', iter, norm(rp), norm(rd), pdgap, tstep);
        
        % 3. update iterates
        x  = x - tstep*dx;
        nu = nu - tstep*dnu;
    end
    
    if iter == MAXITER
        fprintf('Failed to converge in MAXITER=%d iterations\n', MAXITER);
    end
    
    % increase t
    t = t*MU;
end

% evaluate and return objective
fval = 0.5*x'*(qp.P*x) + qp.q'*x;

end


