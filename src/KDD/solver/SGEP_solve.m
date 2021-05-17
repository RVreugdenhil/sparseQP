function [x,F_x] = SGEP_solve(A,B,x0,rho,type,p)

% This function implements the Iteratively Reweighed Quadratic Minorization (IRQM) algorithm (Algorithm 1) described in 
% Song, J.; Babu, P.; Palomar, D.P., "Sparse Generalized Eigenvalue Problem Via Smooth Optimization," IEEE Transactions on Signal Processing, vol.63, no.7, pp.1627-1642, April1, 2015

% It solves the sparse generalized eigenvalue problem via iteratively
% majorizing the penalty function with weighted L2-norm

% Input:
% A: n*n symmetric matrix
% B: n*n symmetric positive definite matrix
% x0: initial point
% rho: regularization parameter
% type: the surrogate function type, including: 'Lp','log' and 'exp'
% p: 0<p<=1 for 'Lp' or p>0 for 'log' and 'exp'

% Output:
% x: the output
% F_x: the objective values at each iteration

MaxIter = 30; % maximum number of iterations
tol = 1e-5; % tolerance, stopping criterion

epsi = 1e-8; % smoothing parameter
n = size(A,1);
norm_diag_A = norm(diag(A));

% Initialize
x = x0;
F_x = zeros(MaxIter,1); % record the objective at each iteration

for k = 1:MaxIter
      
    switch type
        case 'Lp'
            % obj for Lp-norm
            F_x(k) =  x'*A*x - rho*sum(abs(x).^p);
            % weights for L_p norm
            w = p/2*epsi^(p-2)*ones(n,1);
            ind = (abs(x) > epsi);
            w(ind) = (p/2).*(abs(x(ind)).^(p-2));
        case 'log'
            % obj for log
            F_x(k) = x'*A*x - rho/log(1+1/p)*sum(log(abs(x)/p+1));
            % weights for log
            w = (1/(2*epsi*(p+epsi)*log(1+1/p)))*ones(n,1);
            ind = (abs(x) > epsi);
            w(ind) = (0.5/log(1+1/p))./(x(ind).^2+p*abs(x(ind)));
        case 'exp'
            % obj for exp
            F_x(k) = x'*A*x - rho*sum(1-exp(-abs(x)./p));
            % weights for exp
            w = (exp(-epsi/p)/(2*p*epsi))*ones(n,1);
            ind = (abs(x) > epsi);
            w(ind) = exp(-abs(x(ind))./p)./(2*p*abs(x(ind)));
        otherwise
            error('Undefined surrogate function!');
    end

    % Stopping criterion
    if k > 1
        rel_change = abs(F_x(k) - F_x(k-1))/max(1,abs(F_x(k-1))); % relative change in objective
        if rel_change <= tol
            F_x = F_x(1:k);
            return;
        end
    end
    
    % Precondition
    if rho*norm(w)/norm_diag_A > 1e2
        pre_cond = 1./(abs(diag(A))+rho*w);
    else
        pre_cond = ones(n,1);
    end
    
    % Compute the leading generalized eigenvector of the matrix pair (A-rho*diag(w),B)
    
    % steepest ascent method with preconditioning and accerleration (SQUAREM)
    % [x,lambda] = acc_steepest_ascent(A-rho*diag(w),B,x,pre_cond);

    % the preconditioned conjugate gradient method LOBPCG
    % It is usually more stable compared with the preconditioned steepest ascent method



    [x_block,lambda]=lobpcg(x,rho*diag(w)-A,B,@(x)pre_cond.*x,1e-3,10,0);
    x = x_block(:,1);

end