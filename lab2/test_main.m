% test main
% test linear system Ax=b using non-linear cg solver
n = 60;
u = rand(n,1);
Q = RandomRotation(n);
D = diag(u);
A = Q*D*Q';
x = randn(n,1);
b = A*x;

% linear cg solver
tic, x = cg(A,b); toc
norm(A*x-b)

% non-linear cg solver
aux.itr = 1e2;          % iteration limit
aux.tol = 1e-20;         % smallest update of estimate
aux.A   = A;
x0 = zeros(n,1);         % initial guess

fun = @(x)(1/2*(x'*A*x-b'*x));    % loss
grad = @(x)(A*x - b);  % gradient of loss
tic, x_nl = cg_nl(fun,grad,x0,aux); toc
