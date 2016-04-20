function [fval, x, ws, status] = boxqp(qp, varargin)
%BOXQP Solve a sparse quadratic program with factorization caching. 
%   [fval, x, ws, status] = BOXQP(qp, sd, ws0) solves the quadratic program
%
%   minimize   (1/2)x'*P*x + q'*x
%   subject to A*x = b,     (dual var y)
%              G*x + s = h, (dual var z >= 0)
%              s >= 0
%
%   qp: structure containing QP data matrices
%      P, q, A, b, G, h
%   sd: structure containing sparse description
%      p (KKT permutation)
%   ws0: structure containing the warm start point
%      x, s, z, y
%   
%   Other ways to call include
%   ... = BOXQP(qp, sd, ws0) fastest for repeat solutions
%   ... = BOXQP(qp, sd) initial centering
%   ... = BOXQP(qp, [], ws0) no factorization caching
%   ... = BOXQP(qp, [], []) initial centering + no factorization caching
%   ... = BOXQP(qp) is the same as BOXQP(qp, [], [])
%   

% determine calling style
narginchk(1,3);
sd = []; ws0 = [];
if nargin >= 2
    sd = varargin{1};
end
if nargin >= 3
    ws0 = varargin{2};
end

% extract dimensions
nx = size(qp.P, 1); % number of variables
ns = size(qp.G, 1); % number of inequality constraints
nz = ns;
ny = size(qp.A, 1); % number of equality constraints

% algorithm parameters
MU  = 10;           % MU > 1
GAPTOL = 1e-4;      % GAPTOL > 0, require duality gap m/t < GAPTOL
FEASTOL = 1e-4;     % FEASTOL > 0, require |rpri|, |rdual| < FEASTOL
ALPHA = 0.01;       % typically 0.01 to 0.1
BETA  = 0.5;        % typically 0.3 to 0.8
DELTA = 1e-2;       % added to make KKT system sqd
MAXITER = 100;      % maximum number of Newton steps

% determine initial point by solving a closed form problem
if isempty(ws0)
    % check structural invertability of KKT
    % TODO: actually check the rank
    assert(sprank(qp.A) == ny);
    assert(sprank([qp.P, qp.G', qp.A']) == nx);
    
    xzy0 = [ qp.P,         qp.G',         qp.A';
             qp.G, -speye(nz,nz), sparse(nz,ny);
             qp.A, sparse(ny,nz), sparse(ny,ny) ] \ ...
           [ -qp.q;
              qp.h;
              qp.b ];
    x = xzy0(1:nx);
    z = xzy0((nx+1):(nx+nz));
    y = xzy0((nx+nz+1):(nx+nz+ny));
    
    ws0 = struct;
    ws0.x = x;
    ws0.y = y;
    % TODO: determine s0, z0
    s = -z;
    
end
assert(~isempty(ws0));

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
while ny/t > GAPTOL
    % Centering: minimize t*f0(x) + phi(x) 
    %            subject to Ax = b
    fprintf('Centering step\n');
    
    for iter=1:MAXITER
        % 1. compute primal and dual Newton steps
        z  = qp.d - qp.C*x; h = 1./z;
        rp = qp.A*x - qp.b;
        rd = t*(qp.P*x + qp.q) + qp.C'*h + qp.A'*nu;
        pdgap = norm([rp;rd]);
        
        if pdgap <= GAPTOL
            break;
        end
        
        H   = [t*(qp.P) + qp.C'*spdiags(h.^2,0,nz,nz)*qp.C, qp.A'; qp.A, sparse(ny, ny)];
        dH  = [DELTA*speye(nx), sparse(nx,ny); sparse(ny,nx), -DELTA*speye(ny)];
        rhs = -[rd; rp];

        sol = H\rhs;
        %HpdH = H + dH;
        %sol  = ldlsparse(HpdH, sd.p, rhs);
        %sol = sol + ldlsparse(HpdH, sd.p, rhs - H*sol); % refinement step
        dx = sol(1:nx); dnu = sol(nx+1:nx+ny);

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


