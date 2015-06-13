function alpha = analyticlineSearch_cg(A,B,C,x,p)
% compute the analytic solution of line search
% for problem set argmin_x 1/2 || AX - XB - C ||^2_F
% at point x in the search direction p
% the analytic solution of line search is : 
% alpha = - tr((A*x-x*B-C)*(A*p - p*B)^T) / tr((A*p - p*B)*(A*p - p*B)^T)

multi1 = A*x-x*B-C;                     % reg : A*x-x*B-C
multi2 = A*p - p*B;                     % reg : A*p

num = trace(multi1*multi2');            % tr((A*x-x*B-C)*(A*p - p*B)^T)
den = trace(multi2*multi2');            % tr((A*p - p*B)*(A*p - p*B)^T)
alpha = - num / den;


end