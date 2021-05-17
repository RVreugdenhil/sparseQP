function [x,MSE] = beckzero(Q,c,x,eta,sparsity)
% This function implements Algorithm 6 in Beck et al. (write here the
% details of the paper)

    while 1 < 2
        xk = x;
        n = length(c);
        I1 = find(x);
        I0 = 1:n;
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
        [~,Ix] = sort(-p);
        I1(i) = I0(j);
        i = 1;
        while nnz(unique(I1)) < sparsity
            I1(size(I1)+1) = Ix(i);
            i = i+1;
        end
        % step 5 in Algorithm 6
        I1 = unique(I1);

        % step 6 in Algorithm 6
        x = zeros(n,1);
        x(I1) = -0.5*inv(eye(length(I1))/eta+Q(I1,I1))*c(1,I1)';

        if x'*Q*x+c*x + 1/eta*(x'*x) >= xk'*Q*xk+c*xk + 1/eta*(xk'*xk)
            % terminating criterion, line 8 in Algorithm 6
            x = xk;
            MSE = x'*Q*x+c*x + 1/eta*(x'*x);
            break
        end
    end
end


