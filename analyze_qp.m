function sd = analyze_qp(qp)

P = qp.P;
A = qp.A;
C = qp.C;

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
sd = struct;
sd.p     = p{orderingidx};       % best KKT permutation
sd.order = pnames{orderingidx};  % best KKT permutation name
sd.L     = Lfacts{orderingidx};  % sparsity pattern of Cholesky factor
sd.flops = flopcount;

% visualize best structure
figure(1);
subplot(133);
spy(tril(H(sd.p,sd.p)));
hold on;
spy(tril(sd.L ~= 0) ~= tril(H(sd.p,sd.p) ~= 0), 'rx');
hold off;
title(['Cholesky factor with ', sd.order, ', flopcount=', num2str(sd.flops)]);

% for the second two figures, nz denotes fill
end