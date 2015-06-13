function x = cg(A,b,x_start,tol,iter)
% example
%{
% eigenvalues are btw 1 ~ several hundreds(big range)
  n = 600;
  m = 800;
  A = randn(n,m);
  A = A * A';
  b = randn(n,1);
  tic, x = cg(A,b); toc
  norm(A*x-b)
%}
%{
% eigenvalues are btw 0 ~ 1(tight range)
  n = 60;
  u = rand(n,1);
  Q = RandomRotation(n);
  D = diag(u);
  A = Q*D*Q';
  x = randn(n,1);
  b = A*x;
  tic, x = cg(A,b); toc
  norm(A*x-b)
%}

if nargin < 3
    x_start = zeros(size(b));
end
if nargin < 4
    tol = 10^-10;
end
if nargin < 5
    iter = 100;
end

x = x_start;
r = A*x - b;
p = -r;
epsl = 1e-30; % numerical stability

err= [];
for k = 1:numel(b)
    err = [err,norm(r)];
    figure(2), plot(err,'r'), drawnow,hold on, set(gca,'Yscale','log');
    if k == iter + 1
        break
    end
    if norm(r) < tol
        disp('==> solution found!')
        break
    else
        q = A*p; % A*p register
        alpha = ((p'*q) + epsl)\(r'*r);
        x = x + alpha*p;
        r_1 = r;
        r = r_1 + alpha*q;
        beta = ((r_1'*r_1) + epsl)\(r'*r);
        p_1 = p;
        p = -r + beta*p_1;
        orth = p_1'*A*p
    end
    
end

end