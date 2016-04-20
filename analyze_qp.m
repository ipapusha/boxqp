function sd = analyze_qp(qp,varargin)

% check arguments
narginchk(1,2);
quiet = false;
if nargin > 1
    s = varargin{1};
    if strcmpi(s,'quiet')
        quiet = true;
    end
end

P = qp.P;
A = qp.A;
C = qp.C;

% KKT matrix structure
%H = [P + C'*C, A'; A, sparse(size(A,1), size(A,1))];
H = [P + C'*C + 0.1*speye(size(P,1)), A'; A, -0.1*speye(size(A,1))];

if ~quiet
    % visualize KKT sparsity pattern
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

% choose a fill reducing permutation
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

if ~quiet
    % visualize best structure
    % for the second two figures, nz denotes fill
    subplot(133);
    spy(tril(H(sd.p,sd.p)));
    hold on;
    spy(tril(sd.L ~= 0) ~= tril(H(sd.p,sd.p) ~= 0), 'rx');
    hold off;
    title(['Cholesky factor with ', sd.order, ', flopcount=', num2str(sd.flops)]);
end

% TODO: generate initial feasible point by a phase 1 method
% find analytic center of Ax=b and Cx<=d


end





















