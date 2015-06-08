function x = gd(A,x0,aux)
% gradient descent solver
% loss function : argmin_x || AX - I ||^2_F

if ~exist('aux','var'); 
    aux.itr = 100; 
    aux.tol = 1e-6;
end

x = x0;
aux.A = A;

for k = 1 : aux.itr
    
    p = A' - A'*A*x;               % gradient descent direction
    p = p./(norm(p,'fro') + eps); 

    alpha = armijo(x,p,aux);   % step length
    
    chg = alpha*p;

    x = x + chg;                   % next guess
    
    if norm(chg,2) < aux.tol
        sprintf('the number of iterations : %d.', k)
        return
    end
        
end

end