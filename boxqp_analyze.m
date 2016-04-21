function sd = boxqp_analyze(qp,varargin)
%BOXQP_ANALYZE performs offline analysis of KKT factorization matrix
%   sd = BOXQP_ANALYZE(qp, ws0) returns structure description 
%        for use with BOXQP
%
%   qp: structure containing QP data matrices
%      P, q, A, b, G, h
%   ws0: structure containing the warm start point
%      x, s, z, y
%   sd: structure containing sparse description
%      p (KKT permutation)
%
%   Other ways to call include
%   ... = BOXQP_ANALYZE(qp) determines its own starting point ws0
%   ... = BOXQP_ANALYZE(..., 'quiet') does not print or view KKT structure
%
%   See also BOXQP_INITIAL, BOXQP.
% 

% determine calling style
narginchk(1,3);

% determine if we should print progress
quiet = false;
quietidx = strcmp('quiet', varargin);
if any(quietidx)
    quiet = true;
    varargin(quietidx) = [];
end

% figure out parameters
ws0 = [];
if numel(varargin) >= 1
    ws0 = varargin{1};
end

% extract dimensions
nx = size(qp.P, 1); % number of variables
ns = size(qp.G, 1); % number of inequality constraints
nz = ns;
ny = size(qp.A, 1); % number of equality constraints

% algorithm parameters
DELTA = 1e-2;       % DELTA > 0, added to make KKT system sqd

% 1. determine starting point
if isempty(ws0)
    if ~quiet
        ws0 = boxqp_initial(qp);
    else
        ws0 = boxqp_initial(qp,'quiet');
    end;
end
x = ws0.x;
s = ws0.s;
z = ws0.z;
y = ws0.y;

% 2. determine KKT structure
if ~quiet; fprintf('Determining KKT structure.\n'); end;
H = [          qp.P,         sparse(nx,ns),         qp.G',         qp.A';
      sparse(ns,nx), spdiags(z./s,0,ns,ns),  speye(ns,nz), sparse(ns,ny);
               qp.G,          speye(nz,ns), sparse(nz,nz), sparse(nz,ny);
               qp.A,         sparse(ny,ns), sparse(ny,nz), sparse(ny,ny) ];
dH = [  DELTA*speye(nx+ns), sparse(nx+ns,nz+ny); 
       sparse(nz+ny,nx+ns), -DELTA*speye(nz+ny) ];
H = H + dH;

% visualize KKT sparsity pattern
if ~quiet
    subplot(131);
    spy(H);
    title('KKT structure');
end

% calculate unoptimized Cholesky stats
[count, ~, ~, ~, L] = symbfact(H,'sym','lower');
flopcount = sum(count.^2);

if ~quiet
    subplot(132);
    spy(tril(H));
    hold on;
    spy(tril(L ~= 0) ~= tril(H ~= 0), 'rx');
    hold off;
    title(['Cholesky factor, flopcount=', num2str(flopcount)]);
end

% 3. choose a fill reducing permutation
if ~quiet; fprintf('Finding fill-reducing permutation.\n'); end;
p = {amd(H), symrcm(H), colamd(H), colperm(H)};
pnames = {'amd', 'symrcm', 'colamd', 'colperm'};
Lfacts = cell(1, length(p));
flopcounts = zeros(1,length(p));

for i=1:length(p)
    [count, ~, ~, ~, Lfacts{i}] = symbfact(H(p{i},p{i}),'sym','lower');
    flopcounts(i) = sum(count.^2);
end
[flopcount, orderingidx] = min(flopcounts);

% return structure stats
sd = struct;
sd.p     = p{orderingidx};       % best KKT permutation
sd.order = pnames{orderingidx};  % best KKT permutation name
sd.L     = Lfacts{orderingidx};  % sparsity pattern of Cholesky factor
sd.flops = flopcount;

% visualize best structure
% for the second two figures, nz denotes fill
if ~quiet
    subplot(133);
    spy(tril(H(sd.p,sd.p)));
    hold on;
    spy(tril(sd.L ~= 0) ~= tril(H(sd.p,sd.p) ~= 0), 'rx');
    hold off;
    title(['Cholesky factor with ', sd.order, ', flopcount=', num2str(sd.flops)]);
end

end