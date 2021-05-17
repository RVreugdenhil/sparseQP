function [x_br] = determine_M(Q,c,eta, sparsity,stepsize, k, max_iter,period, eigtol)
% Using the dual gradient ascent algorithm to get M
% Input:
% Q, c, eta, sparsity: parameters of original problem
% k: number of eigenvalues used
% eigtol: tolerance for PCA algorithm
% maxIter: maximum iterations

    if nargin < 9
        eigtol = 1e-12; 
    end
    
    %Calculate the eigenvalues and the eigenvalue of matrix Q
    [V,D] = eigs(Q,k,'largestreal','Tolerance',eigtol);
    lambda = diag(D);
    
    
    % get dimensions
    n = length(c);
    
    %Initializing different variables
    alpha = zeros(k, 1);
    Dkbest = -inf;
    for j = 1:max_iter
        
        % stepsize kappa in the paper
        kappa = stepsize/sqrt(j);
        
        % Calculate the subgradient 
        z = zeros(n, 1);
        gamma = (c + alpha'*diag(sqrt(lambda))*V'); 
        
        %take the maximum values of gamma
        [~,idx] = maxk(abs(gamma),sparsity); 
        z(idx) = 1;
        
        % subgradient of alpha
        alphasub = -alpha/2-(eta/2)*diag(sqrt(lambda))*V'*diag(z)*(c + alpha'*diag(sqrt(lambda))*V')'; 
        
        % and perform gradient ascent update
        alpha = alpha + kappa*alphasub;
        
        % take the last p_dp vectors z
        if j>(max_iter-period)
            In(:, j - max_iter + period) = idx;
        end
    end
    x_br = zeros(n, 1);
    x_br(idx) = -0.5*(inv((1/eta)*eye(sparsity)+Q(idx,idx))*c(1,idx)');
    
end