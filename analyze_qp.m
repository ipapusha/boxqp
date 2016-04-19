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

%% prepare MPC

% 1. form objective 1/2*x'*P*x +q'*x
nx = (T+1)*n;
nu = (T)*m;

% i1 = reshape(1:(T*n),n,T);
% i2 = repmat(i1,n,1);
% j1 = repmat(1:(T*n),n,1);
% QQ = sparse(i2(:), j1(:), repmat(Q(:),T,1));
% i1 = reshape(1:(T*m),m,T);
% i2 = repmat(i1,m,1);
% j1 = repmat(1:(T*m),m,1);
% RR = sparse(i2(:), j1(:), repmat(R(:),T,1));

i1 =           reshape(1:(T*n),n,T); % blkdiag(Q,...,Q)
i2 =     T*n + reshape(1:(1*n),n,1); % Qf
i3 = (T+1)*n + reshape(1:(T*m),m,T); % blkdiag(R,...,R)
i4 = repmat(i1,n,1);
i5 = repmat(i2,n,1);
i6 = repmat(i3,m,1);
i = [i4(:); i5(:); i6(:)];

j1 =           repmat(1:(T*n),n,1);
j2 =     T*n + repmat(1:(1*n),n,1);
j3 = (T+1)*n + repmat(1:(T*m),m,1);
j = [j1(:); j2(:); j3(:)];

entries = [repmat(Q(:),T,1); Qf(:); repmat(R(:),T,1)];
P = sparse(i,j,entries, nx+nu, nx+nu);
q = sparse(nx+nu, 1);

% 2. form equality constraints AA*x = b
i1 = 1:n;                      % I initial condition
i2 = n + reshape(1:(T*n),n,T); % blkdiag(A,...,A)
i3 = repmat(i2,n,1);
i4 = n + reshape(1:(T*n),n,T); % blkdiag(B,...,B)
i5 = repmat(i4,m,1);
i6 = n + [1:(T*n)];            % -I on upper blk diagonal
i = [i1(:); i3(:); i5(:); i6(:)];

j1 = 1:n;
j2 =           repmat(1:(T*n),n,1); 
j3 = (T+1)*n + repmat(1:(T*m),n,1);
j4 = n+[1:(T*n)];
j = [j1(:); j2(:); j3(:); j4(:)];

entries = [ones(n,1); repmat(A(:),T,1); repmat(B(:),T,1); repmat(-1,T*n,1)];
AA = sparse(i,j,entries,T*n+n,(T+1)*n+T*m);
b = sparse(1:n, ones(1,n), z, n*(T+1), 1);

% 3. form inequality constraints CC*x <= d
CC = spdiags([-ones(nx+nu,1), ones(nx+nu,1)], [-(nx+nu), 0], 2*(nx+nu), nx+nu);
d = [repmat(xmax,T+1,1); repmat(umax,T,1); repmat(-xmin,T+1,1); repmat(-umin,T,1)];
d = sparse(d);

%% translate to box QP
P = P; q = q;
A = AA; b = b;
C = CC; d = d;
clearvars -except P q A b C d;

%% analyze box QP
% KKT matrix structure
H = [P + C'*C, A'; A, sparse(size(A,1), size(A,1))];

% visualize KKT sparsity pattern
figure(1); 
subplot(131);
spy(H);
title('KKT structure');

% calculate unoptimized Cholesky stats
[count, ~, ~, ~, L] = symbfact(H,'sym','lower');
flopcount = sum(count.^2);
figure(1); 
subplot(132);
spy(tril(H));
hold on;
spy(tril(L ~= 0) ~= tril(H ~= 0), 'rx');
hold off;
title(['Cholesky factor, flopcount=', num2str(flopcount)]);

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
qp = struct;
qp.p     = p{orderingidx};       % best KKT permutation
qp.order = pnames{orderingidx};  % best KKT permutation name
qp.L     = Lfacts{orderingidx};  % sparsity pattern of Cholesky factor
qp.flops = flopcount;

% visualize best structure
figure(1);
subplot(133);
spy(tril(H(qp.p,qp.p)));
hold on;
spy(tril(qp.L ~= 0) ~= tril(H(qp.p,qp.p) ~= 0), 'rx');
hold off;
title(['Cholesky factor with ', qp.order, ', flopcount=', num2str(qp.flops)]);

% for the second two figures, nz denotes fill












