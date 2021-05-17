function [x,obj,In,time] = WarmStart(X,Y,Q,c,eta,sparsity,max_iter,period,bigM)
% Solving the dual function given in Bertsimas & van Parys using a best
% response method
% Input:
% X, Y, eta, sparsity, bigM: parameters of original problem
% max_iter: maximum iterations
% period: solutions used to determine Z
% Output:
% estimator: x*
% obj: x*'Qx* + c'x* + 1/eta x*'x*
% In: size of the reduced problem
% time: computational time to compute Z in seconds



    tic
    % Set initial set of indices
    I = 1:sparsity;

    % Dimensions of the datamatrix X
    N = height(Y);
    n = length(c);
    
    K = zeros(n,1);
    for l = 1:max_iter 
        
        % optimal alpha described in the paper by Bertsimas & van Parys
        alpha = (eye(N)-X(:,I)*(eye(sparsity)/eta+X(:,I)'*X(:,I))^-1*X(:,I)')*Y;
            for i = 1:n
                
            %alpha^T K_i alpha simplified 
            K(i) = norm(alpha'*X(:,i))^2; 
            end
            
        %Select the largest s elements    
        [~,I] = maxk(K,sparsity); 
        
        % save the last p vectors z as desribed in 3.3 in the main paper
            if l>(max_iter-period)
                In(:, l - max_iter + period) = I;
            end
    end
    
    %Only take the unique indices mostly ~ 2 x sparsity
    In = unique(In);  
    time = toc;
    
    % solve P_Z, mentioned in 3.3 in the main paper
    [x_In,obj] = Psolve(Q(In,In),c(1,In),eta,sparsity,bigM);
    x = zeros(length(c),1);
    x(In) = x_In;
end