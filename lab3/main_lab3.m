% lab3
% manifold optimization
% page rank of google

% problem settings
p = 1;
A = [1/5 0 1 0 0;
    1/5 0 0 1/2 0;
    1/5 1/2 0 0 0;
    1/5 0 0 0 1;
    1/5 1/2 0 1/2 0
    ];

n = size(A,1);      % dim of the problem
I = eye(n);
B = 1/n * ones(n);
M = (1-p)*A + p*B;

% solver setting
fun = @(x)(1/2*x'*(M'*M)*x - x'*M'*x + 1/2*x'*x);
grad = @(x)(I - x*x') * (M - I)' * (M - I) * x;
x0 = ones(n,1)./norm(ones(n,1));
% aux.iter = 1e2;
% aux.tol = 1e-6;

% manifold gd solver
x_gd = gd_mani(fun,grad,x0);

% eig, unittest
[V,D] = eigs(M);
x_true = abs(V(:,1))./norm(abs(V(:,1)));

diff = x_true - x_gd