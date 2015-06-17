% lab 2 - non convex opti
%% generate the problem set Ax = b
% quadratic loss function
close all
n = 3;
c = 1;                      % magnitude of eigenvalues

evalue = c*rand(n,1);       % eigenvalue of A
M = diag(evalue);
b = rand(n,1);              % b
Q = RandomRotation(n);

A = Q'*M*Q;                 % pos def sym A
% solver
aux.itr = 1e3;              % iteration limit
aux.tol = 1e-20;            % smallest update of estimate
aux.A   = A;
x0 = zeros(n,1);            % initial guess

fun = @(x)(1/2*(x'*A*x-b'*x));      % loss
grad = @(x)(A*x - b);               % gradient of loss

x_cg = cg(A,b);                         % linear CG solver 
keyboard
x_nlcg = cg_nl(fun,grad,x0,aux);        % non-linear CG solver  

diff_cg = norm(A*x_cg - b,'fro')
diff_nlcg = norm(A*x_nlcg - b,'fro')
keyboard

%% generate the problem set argmin_x 1/2 || AX - XB - C ||^2_F
clear all
A = [1 -2 -2;
    -2 1 -1;
    2 -2 0];
B = [-13 14 -2;
    -10 11 -2;
    5 -4 4];
C = [-5 6 1;
    3 -2 2;
   	-1 -2 -8];

sa = size(A);

% solver - conjugate gradient
% cancellation criteria
aux.itr = 1e1;              % iteration limit
aux.tol = 1e-20;            % smallest update of estimate
aux.A   = A;
aux.B   = B;
aux.C   = C;
% conjugate gradient solver
x0 = zeros(sa);         % initial guess
fun = @(x)(1/2*norm((A*x-x*B-C),'fro')^2);               % loss
const = - A'*C + C*B';
grad = @(x)(A'*A*x + x*B*B' - A'*x*B - A*x*B' + const);  % gradient of loss

x_cg = cg_nl(fun,grad,x0,aux);          % conjugate gradient solver       
% x_gd = gd(fun,grad,x0,aux);             % gradient descent solver

% check sol 
% diff = norm(A*x_cg-x_cg*B-C,'fro')
% diff = norm(A*x_gd-x_gd*B-C,'fro')
