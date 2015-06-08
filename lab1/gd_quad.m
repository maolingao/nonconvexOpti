function x = gd_quad(A,b,x0,aux)
% gradient descent solver
% for quadratic function f(x) = 1/2x^TAx - b^Tx

if ~exist('aux','var'); 
    aux.itr = 100; 
    aux.tol = 1e-10;
    aux.tolgrad = 1e-6;
end

x = x0;
 
for k = 1 : aux.itr
    
    d = A*x - b;                   % gradient descent direction

    t = (d'*d) / (d'*A*d + eps);    % step length
    
    chg = t*d;

    x = x - chg;                   % next guess
    
    if norm(chg,2) < aux.tol || norm(d,2) < aux.tolgrad
        sprintf('the number of iterations : %d.', k)
        disp('the minimum :')
        x
        return
    end
        
end

end