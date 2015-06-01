% ex 1 - non convex opti
%% generate the problem set Ax = b
n = 2;
c = 1; % magnitude

evalue = c*rand(n,1); % eigenvalue of A
M = diag(evalue);
b = rand(n,1); % b

Q = RandomRotation(n);

A = Q'*M*Q;     % matrix A

%% solver - gradient descent

% cancellation criteria
aux.itr = 1e3;
aux.tol = 1e-11;
aux.tolgrad = 1e-9;

% gradient descent
x0 = zeros(n,1);        % initial guess
x_gd = gd(A,b,x0,aux);  

% true solution
x_true = A\b;           

% check sol 
diff = max(abs(x_gd - x_true))
