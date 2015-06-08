% ex 1 - non convex opti
%% generate the problem set argmin_x || AX - I ||^2_F

A = [2 3 0 1;
    0 7 3 0;
    1 3 1 1;
    1 1 0 1];
sa = size(A);
%% solver - gradient descent

% cancellation criteria
aux.itr = 1e5;
aux.tol = 1e-6;

% gradient descent
x0 = zeros(sa);        % initial guess
x_gd = gd(A,x0,aux);   % gradient descent solver

% true solution
x_true = A\eye(sa);           

% check sol 
diff = max(max((abs(x_gd - x_true))))
