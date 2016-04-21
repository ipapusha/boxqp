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
%   ... = BOXQP(..., 'quiet') does not print anything during execution
%
%   See also BOXQP_INITIAL, BOXQP_ANALYZE.
%   

% determine calling style
narginchk(1,4);

% determine if we should print progress
quiet = false;
quietidx = strcmp('quiet', varargin);
if any(quietidx)
    quiet = true;
    varargin(quietidx) = [];
end

% figure out parameters
sd = []; ws0 = [];
if numel(varargin) >= 1
    sd = varargin{1};
end
if numel(varargin) >= 2
    ws0 = varargin{2};
end

% extract dimensions
nx = size(qp.P, 1); % number of variables
ns = size(qp.G, 1); % number of inequality constraints
nz = ns;
ny = size(qp.A, 1); % number of equality constraints

% algorithm parameters
MU  = 10;           % MU > 1
GAPTOL = 1e-5;      % GAPTOL > 0, require duality gap m/t < GAPTOL
RGAPTOL = 1e-10;    % RGAPTOL > 0, used to check infeasibility
FEASTOL = 1e-5;     % FEASTOL > 0, require |rpri|, |rdual| < FEASTOL
ALPHA = 0.01;       % 0 < ALPHA < 1/2, typically 0.01 to 0.1
BETA  = 0.5;        % 0 < BETA < 1, typically 0.3 to 0.8
DELTA = 1e-2;       % DELTA > 0, added to make KKT system sqd
MAXITER = 100;      % maximum number of Newton steps

% determine initial if necessary
if isempty(ws0)
    if ~quiet
        ws0 = boxqp_initial(qp);
    else
        ws0 = boxqp_initial(qp,'quiet');
    end;
end

% check strict feasibility of initial point
assert(~isempty(ws0));
assert(~isempty(ws0.x));
assert(~isempty(ws0.y));
assert(all(ws0.s >= 0));
assert(all(ws0.z >= 0));

if ~quiet; 
    fprintf('===============================================\n'); 
    fprintf('=== BOXQP primal dual path following method ===\n'); 
    fprintf('===============================================\n'); 
end;

% primal dual path following method
x = ws0.x;
s = ws0.s;
z = ws0.z;
y = ws0.y;
for iter=1:MAXITER
    % 1. Evaluate residuals, gap, and stopping criteria
    gap = z'*s;
    mu = gap/(MU*nz);  %t^(-1) = (MU*nz/gap)^(-1);
    rx = qp.P*x + qp.q + qp.G'*z + qp.A'*y;
    rs = (z.*s) - mu;
    rz = qp.G*x + s - qp.h;
    ry = qp.A*x - qp.b;
    
    % quit if tolerances met
    if ((norm([rz; ry], 2) <= FEASTOL) && ... % primal feas
        (norm(rx, 2) <= FEASTOL) && ...       % dual feas
        (gap <= GAPTOL))                      % primal - dual
        break;
    end
    
    % 2. Compute primal-dual search direction
    H = [          qp.P,         sparse(nx,ns),         qp.G',         qp.A';
          sparse(ns,nx), spdiags(z./s,0,ns,ns),  speye(ns,nz), sparse(ns,ny);
                   qp.G,          speye(nz,ns), sparse(nz,nz), sparse(nz,ny);
                   qp.A,         sparse(ny,ns), sparse(ny,nz), sparse(ny,ny) ];
    rhs = [-rx; mu./s - z; -rz; -ry];
    
    if isempty(sd)
        % exact search direction
        dxszy = H\rhs;
    else
        % approximate search direction with 
        % factorization and iterative refinement
        dH = [  DELTA*speye(nx+ns), sparse(nx+ns,nz+ny); 
               sparse(nz+ny,nx+ns), -DELTA*speye(nz+ny) ];
        HpdH = H + dH;
        rhs = full(rhs);
        dxszy = ldlsparse(HpdH, sd.p, rhs);
        dxszy = dxszy + ldlsparse(HpdH, sd.p, rhs - H*dxszy); % refinement step
    end
    
    dx = dxszy(1:nx);
    ds = dxszy((nx+1):(nx+ns));
    dz = dxszy((nx+ns+1):(nx+ns+nz));
    dy = dxszy((nx+ns+nz+1):(nx+ns+nz+ny));
    
    % 3. Modified backtracking line search
    s1 = maxstep(z, dz);
    s2 = maxstep(s, ds);
    step = 0.99*min(s1, s2);
    
    % multiply step by BETA until residual is decreased
    rcur = norm([rx; rs; rz; ry], 2);
    while 1
        % new point
        xn = x+step*dx; 
        sn = s+step*ds; 
        zn = z+step*dz; 
        yn = y+step*dy;
        
        % new gap
        gapn = zn'*sn;
        mun = gapn/(MU*nz);
        
        % new residuals
        rxn = qp.P*xn + qp.q + qp.G'*zn + qp.A'*yn;
        rsn = (zn.*sn) - mun;
        rzn = qp.G*xn + sn - qp.h;
        ryn = qp.A*xn - qp.b;
        
        % quit line search if residual is decreased
        if norm([rxn; rsn; rzn; ryn], 2) <= (1 - ALPHA*step)*rcur
            break;
        end
        
        % otherwise, decrease step
        step = BETA*step;
    end
    
    % detect infeasibility
    rgap = abs((gapn - gap)/gap);
    if (rgap <= RGAPTOL) && (gap > GAPTOL)
        break;
    end
    
    % 4. Update iterates
    if ~quiet
        fprintf('%-2d gap: %-11g, rgap: %6.2e, rx: %6.2e, rs: %6.2e, rz: %6.2e, ry: %6.2e\n', ...
                iter, gap, rgap, norm(rx), norm(rs), norm(rz), norm(ry));
    end
    
    x = x+step*dx;
    s = s+step*ds;
    z = z+step*dz;
    y = y+step*dy;
end

% return answers
fval = 0.5*x'*(qp.P*x) + qp.q'*x;
ws = struct;
ws.x = x;
ws.s = s;
ws.z = z;
ws.y = y;

status = struct;
status.iter = iter;
status.desc = 'solved';
if (rgap <= RGAPTOL) && (gap > GAPTOL)
    status.desc = 'infeasible';
elseif (iter == MAXITER)
    status.desc = 'max iterations reached';
end

if ~quiet
    if strcmpi(status.desc, 'infeasible')
        fprintf('relative gap decrease: %g, suspect infeasibility\n', rgap);
    elseif strcmpi(status.desc, 'max iterations reached')
        fprintf('Failed to converge in MAXITER=%d iterations\n', MAXITER);
    end
end
status.gap = gap;
status.rx = norm(rx);
status.rs = norm(rs);
status.rz = norm(rz);
status.ry = norm(ry);

end



function step = maxstep(w, dw)
% step = MAXSTEP(w, dw)
% returns sup(t in [0,1] | w + step*dw >= 0)
%         = min(1, min_i(-x_i/dx_i | dx_i < 0))
if all(dw >= 0)
    step = 1;
else
    neg = (dw < 0);
    step = full(min(1, min(-w(neg)./dw(neg))));
end

end

