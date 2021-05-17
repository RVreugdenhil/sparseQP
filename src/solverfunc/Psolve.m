function [x_P, obj, I_P] = Psolve(Q,c,eta,sparsity, bigM, ops)
% Solve problem (P_Z) defined in Section 3.3 in the main paper
% Input:
% Q, c, eta, sparsity, bigM: parameters of original problem
% ops: options that can be specified for the MOSEK solver
% Output:
% estimator: x*
% obj: Jk* defined in the paper
% I_P: indices corresponding to the reduced problem

    if nargin < 6
        ops = sdpsettings('solver','mosek','mosek.MSK_DPAR_OPTIMIZER_MAX_TIME', 300,'verbose',0);
    end
    
    yalmip('clear');
    
    % get dimensions
    n = length(c);
    
    % solve P_Z using an MIQP solver
    x = sdpvar(n,1);
    z = binvar(n,1);
    residuals = x'*Q*x+c*x+(1/eta)*(x'*x);
    Constraint = [sum(z) <= sparsity,-bigM*z <= x, x <= bigM*z];
    
    output = optimize(Constraint,residuals,ops);
    
    x_P = value(x);
    I_P = find(value(z));
    obj = value(residuals);
end