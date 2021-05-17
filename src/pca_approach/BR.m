function [x,obj,In,time] = BR(Q,c,eta, sparsity,k, max_iter,period, bigM,eigtol)
% Best response algorithm
% Input:
% Q, c, eta, sparsity, bigM: parameters of original problem
% k: number of eigenvalues used
% eigtol: tolerance for PCA algorithm
% maxIter: maximum iterations
% period: solutions used to determine Z
% Output:
% estimator: x*
% obj: x*'Qx* + c'x* + 1/eta x*'x*
% In: size of the reduced problem
% time: computational time to compute Z in seconds

    if nargin < 9
        eigtol = 1e-12; 
    end
    
    tic 
    
    %Calculate the eigenvalues and the eigenvalue of matrix Q
    [V,D] = eigs(Q,k,'largestreal','Tolerance',eigtol);
    lambda = diag(D);

    % get dimensions
    n = length(c);

    %Initializing z
    z = zeros(n,1);
    for j = 1:max_iter
        
        % Best response on alpha_t
        alpha = -inv(eye(k)/eta+diag(sqrt(lambda))*V'*diag(z)*V*diag(sqrt(lambda)))*diag(sqrt(lambda))*V'*diag(z)*c';
        
        % determine z_t 
        gamma = c + alpha'*diag(sqrt(lambda))*V'; 
        %take the s largest absolute values of gamma
        [~,idx] = maxk(abs(gamma),sparsity); 
        z = zeros(n, 1);
        z(idx) = 1;
        
        % take the last p_br vectors z
        if j>(max_iter-period)
            In(:, j - max_iter + period) = idx;
        end
    end
    
    % make vector Z, ||Z||_0 ~2 x sparsity
    In = unique(In); 
    time = toc;
    
    %solve P_Z mentioned in Section 3.3
    [x_In,obj] = Psolve(Q(In,In),c(1,In),eta,sparsity,bigM);
    x = zeros(size(Q,1),1);
    x(In) = x_In;
    
end