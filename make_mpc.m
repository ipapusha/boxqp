function qp = make_mpc(A, B, z, xmin, xmax, umin, umax, Q, Qf, R, T)

% extract dimensions
n = size(A,1);
m = size(B,2);

% set up the QP description structure
qp = struct;

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
qp.P = sparse(i,j,entries, nx+nu, nx+nu);
qp.q = sparse(nx+nu, 1);

% 2. form equality constraints A*x = b
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
qp.A = sparse(i,j,entries,T*n+n,(T+1)*n+T*m);
qp.b = sparse(1:n, ones(1,n), z, n*(T+1), 1);

% 3. form inequality constraints CC*x <= d
qp.C = spdiags([-ones(nx+nu,1), ones(nx+nu,1)], [-(nx+nu), 0], 2*(nx+nu), nx+nu);
qp.d = sparse([repmat(xmax,T+1,1); repmat(umax,T,1); repmat(-xmin,T+1,1); repmat(-umin,T,1)]);

end