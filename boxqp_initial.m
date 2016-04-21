function ws0 = boxqp_initial(qp,varargin)
%BOXQP_INITIAL calculates an initial point by solving a centering problem
%   ws0 = BOXQP_INITIAL(qp) returns the (analytical) solution to
%
%   minimize   (1/2)x'*P*x + q'*x + (1/2)s'*s
%   subject to A*x = b,     (dual var y)
%              G*x + s = h, (dual var z)
%
%   qp: structure containing QP data matrices
%      P, q, A, b, G, h
%   ws0: structure containing the warm start point
%      x, s, z, y
%
%   Other ways to call include
%   ... = BOXQP_INITIAL(..., 'quiet') does not print anything 
%         during execution
%
%   See also BOXQP_ANALYZE, BOXQP.
% 

% determine calling style
narginchk(1,2);

% determine if we should print progress
quiet = false;
quietidx = strcmp('quiet', varargin);
if any(quietidx)
    quiet = true;
    varargin(quietidx) = [];
end

% extract dimensions
nx = size(qp.P, 1); % number of variables
ns = size(qp.G, 1); % number of inequality constraints
nz = ns;
ny = size(qp.A, 1); % number of equality constraints

% 1. figure out starting point
if ~quiet; fprintf('Calculating initial point.\n'); end;

% check structural invertability of KKT
% TODO: actually check the rank
assert(sprank(qp.A) == ny);
assert(sprank([qp.P, qp.G', qp.A']) == nx);

% minimize   (1/2)x'*P*x + q'*x + (1/2)s'*s
% subject to A*x = b,     (dual var y)
%            G*x + s = h, (dual var z)

xzy0 = [ qp.P,         qp.G',         qp.A';
         qp.G, -speye(nz,nz), sparse(nz,ny);
         qp.A, sparse(ny,nz), sparse(ny,ny) ] \ ...
       [ -qp.q;
          qp.h;
          qp.b ];
x = xzy0(1:nx);
z = xzy0((nx+1):(nx+nz));
s = -z;
y = xzy0((nx+nz+1):(nx+nz+ny));

ws0 = struct;
ws0.x = x;
ws0.y = y;

ap = -min(-z);
if ap < 0
    ws0.s = -z;
else
    ws0.s = -z + (1+ap);
end

ad = -min(z);
if ad < 0
    ws0.z = z;
else
    ws0.z = z + (1+ad);
end

% 2. check strict feasibility of initial point
assert(~isempty(ws0));
assert(~isempty(ws0.x));
assert(~isempty(ws0.y));
assert(all(ws0.s >= 0));
assert(all(ws0.z >= 0));