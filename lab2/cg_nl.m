function x = cg_nl(fun,grad,x0,aux)
% conjugate gradient solver
% loss function : argmin_x 1/2 || AX - XB - C ||^2_F
% fun : loss function handle
% grad : loss function gradient handle
% x0 : initial guess
% aux : cancelation parameters

if ~exist('aux','var'); 
    aux.itr = 100; 
    aux.tol = 1e-6;
end
keyboard
x = x0;
r = grad(x);                          % initial gradient
p = -r;        % initial gradient descent 

for k = 1 : aux.itr

%         q = aux.A*p; % A*p register
%         alpha = ((p'*q) + eps)\(r'*r);
    alpha = armijo(fun,x,r,p);          % step length - armijo
%     alpha = analyticlineSearch_cg(aux.A,aux.B,aux.C,x,p); % analytical line search
    
    chg = alpha*p;

    x = x + chg;                        % next guess
    fun(x) - fun(x-chg)
%     keyboard
%     assert(fun(x) - fun(x-chg) <= 0, 'ATTENTION : loss funtion is increasing!')
    
%     if norm(chg,2) < aux.tol
%         sprintf('the number of iterations : %d.', k)
%         return
%     end
    r_1 = r; 
    r = grad(x);                        % residual
    
    beta = (norm(r_1,'fro')^2+eps) \ norm(r,'fro')^2;             % beta Fletcher-Reeves update rule
    p_1 = p;                            % cash the old p_1
    
    p = -r + beta*p_1;                   % conjugate gradient direction
%     orth = p_1'*aux.A*p
        
end

end