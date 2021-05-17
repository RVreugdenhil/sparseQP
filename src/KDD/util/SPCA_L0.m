function [x,F_x] = SPCA_L0(A,b,x0,rho)

% This function implements Algorithm 4 described in 
% Song, J.; Babu, P.; Palomar, D.P., "Sparse Generalized Eigenvalue Problem Via Smooth Optimization," IEEE Transactions on Signal Processing, vol.63, no.7, pp.1627-1642, April1, 2015

% It solves the following problem:
% max  x'*A*x - rho*||x||_0
% s.t. x'*Diag(b)*x = 1

% Input:
% A: n*n symmetric PSD matrix
% b: n*1 positive vector, b > 0
% x0: initial point
% rho: regularization parameter, rho >= 0

% Output:
% x: the output
% F_x: the objective values at each iteration

MaxIter = 1000;
tol = 1e-5;
n = size(A,1);

% compute the auxiliary matrix A_tilde
b_sqrtinv = 1./sqrt(b);
A_tilde = diag(sparse(b_sqrtinv))*A*diag(sparse(b_sqrtinv));

% Initialize
x = x0;
F_x = zeros(MaxIter,1);

for k = 1:MaxIter
    
     % record the objective
     a = 2*(A_tilde*x);
     F_x(k) = 0.5*(x'*a) - rho*sum(abs(x)>1e-10);    
     
      % Stopping criterion
      if k > 1
        rel_change = abs(F_x(k) - F_x(k-1))/max(1,abs(F_x(k-1))); % relative change in objective
        if rel_change <= tol
            x = b_sqrtinv.*x;
            F_x = F_x(1:k);
            return;
        end
      end
      
    % update x
    x = x_update(a,rho);
end
x = b_sqrtinv.*x;

function x = x_update(a,rho)
% update x according to Proposition 6
x = sparse(length(a),1);
ind_nonzero = find(abs(a) >= rho); % only consider entries with absolute value >= rho 
if isempty(ind_nonzero) 
    [a_val,ind_s] = max(abs(a));
    if a_val ~= 0
        x(ind_s) = sign(a_val);
    else
        x(ind_s) = 1; % for the case a is a vector of all zeros
    end
else
    % find the smallest integer s such that the L2 norm of the largest s+1
    % elements of 'a' minus the L2 norm of the largest s elements of 'a' is smaller
    % or euqal to rho, and return the indices of the s largest (in absolute value) entries
    [b,ind] = sort(a(ind_nonzero).^2,'descend');
    b_norms = sqrt(cumsum(b));
    b_diff = b_norms(2:end) - b_norms(1:end-1);
    s = find([b_diff;0] <= rho,1); % append one 0 at the end to ensure we can find the correct index in extreme cases
    ind_s = ind_nonzero(ind(1:s));
    x(ind_s) = a(ind_s)/norm(a(ind_s));
end