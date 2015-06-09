% ex 2 - non convex opti
%% generate the problem set argmin_x || AX - I ||^2_F
A = [2 3 0 1;
    0 7 3 0;
    1 3 1 1;
    1 1 0 1];
sa = size(A);

%% solver - gradient descent
% cancellation criteria
aux.itr = 1e5;          % iteration limit
aux.tol = 1e-6;         % smallest update of estimate

% gradient descent solver
x0 = zeros(sa);         % initial guess
fun = @(x)(1/2*norm((A*x-eye(sa)),'fro')^2);    % loss
grad = @(x)(A'*A*x - A');                       % gradient of loss

x_gd = gd(fun,grad,x0,aux);                     % gradient descent solver

% true solution
x_true = A\eye(sa);           

% check sol 
diff = norm(x_gd - x_true)
