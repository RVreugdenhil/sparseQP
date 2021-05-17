function [x,obj,In,time] = DP(Q,c,eta, sparsity,k,stepsize,max_iter, period, bigM,eigtol)
% Dual gradient ascent
% Input:
% Q, c, eta, sparsity, stepsize, bigM: parameters of original problem
% k: number of eigenvalues used
% eigtol: tolerance for PCA algorithm
% maxIter: maximum iterations
% period: solutions used to determine Z
% Output:
% estimator: x*
% obj: x*'Qx* + c'x* + 1/eta x*'x*
% In: size of the reduced problem
% time: computational time to compute Z in seconds


    if nargin < 10
        eigtol = 1e-12; 
    end
    
    tic
    [V,D] = eigs(Q,k,'largestreal','Tolerance',eigtol);%Calculate the eigenvalues and the eigenvalue of matrix Q
    lambda = diag(D);
        
    % get dimensions
    n = length(c);
    
    % Initializing alpha
    alpha = zeros(k, 1);
    for j = 1:max_iter
        
        % stepsize kappa_t in the paper
        kappa = stepsize/sqrt(j);
        
        % determine z_t 
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
    
    % make vector Z, ||Z||_0 ~sparsity
    In = unique(In);
    In(In == 0 )= [];
    time = toc;
    
    % solving P_Z as mentioned in Section 3.3 in the main paper
    [x_In,obj] = Psolve(Q(In,In),c(1,In),eta,sparsity,bigM);
    x = zeros(length(c),1);
    x(In) = x_In;
end