function alpha = analyticlineSearch(A,x,p)
% compute the analytic solution of line search
% for problem set argmin_x || AX - I ||^2_F
% at point x in the search direction p
% the analytic solution of line search is : 
% alpha = - tr((A*x-I)(A*p)^T) / tr((A*p)(A*p)^T)

AX_I = A*x - eye(size(A));      % reg : A*x-I
AP = A*p;                       % reg : A*p
num = trace(AX_I*AP');          % tr((A*x-I)(A*p)^T)
den = trace(AP*AP');            % tr((A*p)(A*p)^T)
alpha = - num / den;


end