function [x_lambda] = tikhonovFast(U,s,V,b,lambda)
%TIKHONOV Tikhonov regularization.
%
% [x_lambda,rho,eta] = tikhonov(U,s,V,b,lambda,x_0)
% [x_lambda,rho,eta] = tikhonov(U,sm,X,b,lambda,x_0) ,  sm = [sigma,mu]
%
% Computes the Tikhonov regularized solution x_lambda, given the SVD or
% GSVD as computed via csvd or cgsvd, respectively.  If the SVD is used,
% i.e. if U, s, and V are specified, then standard-form regularization
% is applied:
%    min { || A x - b ||^2 + lambda^2 || x - x_0 ||^2 } .
% If, on the other hand, the GSVD is used, i.e. if U, sm, and X are
% specified, then general-form regularization is applied:
%    min { || A x - b ||^2 + lambda^2 || L (x - x_0) ||^2 } .
%
% If an initial estimate x_0 is not specified, then x_0 = 0 is used.
%
% Note that x_0 cannot be used if A is underdetermined and L ~= I.
%
% If lambda is a vector, then x_lambda is a matrix such that
%    x_lambda = [ x_lambda(1), x_lambda(2), ... ] .
%
% The solution norm (standard-form case) or seminorm (general-form
% case) and the residual norm are returned in eta and rho.

% Per Christian Hansen, DTU Compute, April 14, 2003.

% Reference: A. N. Tikhonov & V. Y. Arsenin, "Solutions of Ill-Posed
% Problems", Wiley, 1977.

% Special case of tikhonov function for cannula microscopy

    if (min(lambda)<0)
      error('Illegal regularization parameter lambda')
    end
    
    n = size(V,1);    
    beta = U'*b;
    zeta = s(:,1).*beta;
    ll = length(lambda);    
    x_lambda = zeros(n,ll);
    
    % Treat each lambda separately.
    for i=1:ll
        x_lambda(:,i) = V*(zeta./(s.^2 + lambda(i)^2));
    end
    
end