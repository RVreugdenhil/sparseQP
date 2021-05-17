function [x,MSE,Iopt] = beckfull(Q,c,x,eta,sparsity)
% This function implements Algorithm 7 in Beck et al. (write here the
% details of the paper)

% step 1 of algorithm 7 by Beck et al. 
% this is similar to the beckzero function which is based on Algorithm 6 of
% the paper
    while 1 < 2
    xk = x;
    n = length(c);
    I1 = find(x);
    I0 = [1:n]';
    I0(I1) = [];
    % calculate p of the negative gradient of f(x^k)
    p = abs(-((eye(length(n))/eta+Q)*x + c'));
    % step 2 in Algorithm 6
    [~,i] = min(p(I1));
    % step 3 in Algorithm 6
    [~,j] = min(-p(I0));
    if isempty(i) == 1
        i = 1;
    end
    % step 4 in Algorithm 6
    I1(i) = I0(j);
    [~,Ix] = sort(-p);
    i = 1;
    while nnz(unique(I1)) < sparsity
        I1(size(I1)+1) = Ix(i);
        i = i+1;
    end
    % step 5 in Algorithm 6
    I1 = unique(I1);
    I0 = [1:n]';
    I0(I1) = [];

    Ior = I1;
    optMSE = inf;
    % step 2,3,4
    for i = 1:size(I1)
        for j = 1:size(I0)
            I1 = Ior;
            I1(i) = I0(j);
            x = zeros(n,1);
            x(I1) = -0.5*inv(eye(length(I1))/eta+Q(I1,I1))*c(1,I1)';
            if x(I1)'*Q(I1,I1)*x(I1) + c(I1)*x(I1) + 1/eta*(x(I1)'*x(I1)) < optMSE
                optMSE = x'*Q*x+c*x + 1/eta*(x'*x);
                Iopt = I1;
            end
        end
    end
    % step 6 in Algorithm 7
    x = zeros(n,1);
    x(Iopt) = -0.5*inv(eye(length(Iopt))/eta+Q(Iopt,Iopt))*c(1,Iopt)';

    if x'*Q*x+c*x + 1/eta*(x'*x) >= xk'*Q*xk+c*xk + 1/eta*(xk'*xk)
        % terminating criterion, line 7 in Algorithm 6
        x = xk;
        MSE = x'*Q*x+c*x + 1/eta*(x'*x);
        break
    end
end
end


