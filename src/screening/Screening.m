function [x,obj,time,sizeIn] = Screening(W,Y,Q,c,eta,sparsity,maxIter,bigM)
% Safe screening method by Atamturk & Gomez using alternating optimization
% over x and z 
% Input:
% W, y, eta, sparsity: parameters of original problem
% k: number of eigenvalues used
% maxIter: maximum iterations
% Output:
% estimator: x*
% obj: x*'Qx* + c'x* + 1/eta x*'x*
% time: computational time to compute Z in seconds
% sizeIn : the reduced problem size and the amount of variables we
% can fix to nonzero


    loop = tic;
    
    % dimension of the data matrix
    N = height(Y);
    n = length(c);
    
    % give initial value to x*
    if n > N
            x_s = (eta*eye(n) - eta*W'*(eye(N)+W*eta*W')^(-1)*W*eta)*W'*Y;
        else
            x_s = (inv((1/eta)*eye(n)+W'*W)*W'*Y);
    end
    
    % optimization of variables z and x which are intercorrelated           
    for i = 1:maxIter
        
        % optimize over z
        z = sdpvar(n,1);
        residuals = x_s'*diag(1./z)*x_s;
        Constraint = [sum(z) <= sparsity, 0 <= z <= 1 ];
        ops = sdpsettings('solver','mosek','verbose',0);
        optimize(Constraint,residuals,ops);
        z_s = value(z);    
        
        % optimize over x
        if n > N
            x_s = (eta*diag(z_s) - eta*diag(z_s)*W'*(eye(N)+W*eta*diag(z_s)*W')^(-1)*W*eta*diag(z_s))*W'*Y;
        else
            x_s = (((1/eta)*diag(1./z_s)+W'*W)^-1*W'*Y);
        end
      
    end

    % solution to the relaxed problem
    AnsScreen = ((Y-W*x_s)'*(Y-W*x_s)+(1/eta)*sum(x_s.^2./(z_s)))/N; 
    
    % determine delta mentioned in the paper by Atamturk & Gomez
    eps = Y-W*x_s;
    delta = (W'*eps).^2;
    [~,I] = maxk(delta,sparsity+1);
    deltak = delta(I(end-1));
    deltak1 = delta(I(end)); 
    
    % round the values of the relaxed solution to get an upperbound to the
    % relaxed problem
    x_relax = -0.5*(inv((1/eta)*eye(sparsity)+Q(I(1:end-1),I(1:end-1)))*c(1,I(1:end-1))');
    upp = ((Y-W(:,I(1:end-1))*x_relax)'*(Y-W(:,I(1:end-1))*x_relax)+(1/eta)*(x_relax'*x_relax))/N;

    % fix the variabels of z to be either 0 or 1, if we do not fix it set it to
    % 99 details are listed in the paper by Atamturk & Gomez
    zn = zeros(n,1);
    for i = 1:n
        if delta(i) <= deltak1 && upp<=AnsScreen - eta*(delta(i)-deltak)
            zn(i) = 0;
        elseif -delta(i) <= -deltak && upp<= AnsScreen + eta*(delta(i)-deltak1)
            zn(i) = 1;
        else
            zn(i) = 99;
        end
    end 
    Iopt = zeros(sparsity,1);
    
    % Optimizing over only the non fixed indices

    % if we fix all variables P_Z is a convex quadratic programming problem
    if size(find(zn==1),1) ==sparsity
        Iopt(1:size(find(zn==1),1)) = find(zn==1);
        sizeIs = 0;
        sizeRed = 0;
        sizeIn = sizeIs + sparsity-sizeRed;
        time = toc(loop);
        x = zeros(size(Q,1),1);
        x(Iopt) = -0.5*(inv((1/eta)*eye(sparsity)+Q(Iopt,Iopt))*c(1,Iopt)');
        obj = x'*Q*x+c*x+(1/eta)*x'*x;
        
    % if we fix at least one variable z_i to be equal to 1
    elseif size(find(zn==1),1) > 0
        Iopt(1:size(find(zn==1),1)) = find(zn==1);
        In = find(zn==99);
        sizeIs = size(find(zn==99),1);
        sizeRed = sparsity - size(find(zn==1),1);
        sizeIn = sizeIs + sparsity-sizeRed;
        time = toc(loop);
        
        %solve P_Z as mention in Section 3.3 in the main paper
        [~,~,I_P] = Psolve(Q(In,In),c(1,In),eta,sizeRed,bigM);
        Iopt(size(find(zn==1),1)+1:end,1) = In(I_P);
        x = zeros(size(Q,1),1);
        x(Iopt) = -0.5*(inv((1/eta)*eye(sparsity)+Q(Iopt,Iopt))*c(1,Iopt)');
        obj = x'*Q*x+c*x+(1/eta)*x'*x;
        
    % no variable z_i is fixed to 1
    else
        sizeIs = size(find(zn==99),1);
        sizeRed = sparsity;
        sizeIn = sizeIs + sparsity-sizeRed;
        In = find(zn);
        time = toc(loop);
        
        %solve P_Z as mention in Section 3.3 in the main paper
        [x_In,obj] = Psolve(Q(In,In),c(1,In),eta,sparsity,bigM);
        x = zeros(length(c),1);
        x(In) = x_In;
    end
end