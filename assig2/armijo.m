function alpha = armijo(x,p,aux,par)
% line-search : armijo condition (backtracking)
if ~exist('par','var')
    par.c = 0.9;
    par.rho = 0.5;
    par.alpha0 = 1;
end

% paremeter to be tuned
c = par.c;
rho = par.rho;
alpha0 = par.alpha0;

% aux information
A = aux.A;
sa  = size(A);
I   = eye(sa);

% initial settings
alpha   = alpha0; 
diff    = 1;
f_x     = 1/2 * norm(A*x-I,'fro')^2;                    % loss function value @ x

while diff > 0
    
    f_func  = 1/2 * norm(A*(x+alpha*p)-I,'fro')^2;      % loss function value
    f_cord  = f_x + alpha * c * norm(p,'fro')^2;        % cord line value
    diff    = f_func - f_cord;                          % compare 
    if diff > 0
        alpha   = alpha * rho;
    end
end

end