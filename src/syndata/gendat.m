function [X,Y,wtrue,Itrue] = gendat(N,n,seed,rho,sparsity,SNR)
% Generate synthetic dataset as described in the paper
% Input:
% N, n, rho, sparsity, SNR: parameters of original problem
% seed used for reproducability
% Output:
% X, Y: Input, Output matrix
% wtrue, Itrue: true solution used to generate the data


    if sparsity >= n
        msg = 'Error! sparsity level is the same or higher than the feature dimension';
        disp(msg)
    else
        rng(seed,'v5uniform');
        
        % create Sigma matrix
        sigma = zeros(n,n);
        for i = 1:n
            for j = 1:n
                sigma(i,j) = rho^(abs(i-j)); 
            end
        end
    
        % generate data matrix X 
        X = mvnrnd(zeros(N,n),sigma,N); 
        
        % index of nonzero components
        Itrue = sort(randperm(n,sparsity))';
        
        % create true vector
        wtrue = zeros(n,1);
        R = [-1,1];
        for h = 1:sparsity
            wtrue(Itrue(h)) = R(randi([1 2]));
        end
        noise = var(X*wtrue)/(SNR);
        Y = X*wtrue+normrnd(0,sqrt(noise),[N,1]); %Output 
    end
end