function [x,F_x] = SPCA_approx(A,b,x0,rho,type,p)

% This function implements Algorithm 3 described in 
% Song, J.; Babu, P.; Palomar, D.P., "Sparse Generalized Eigenvalue Problem Via Smooth Optimization," IEEE Transactions on Signal Processing, vol.63, no.7, pp.1627-1642, April1, 2015

% It solves the following problem:
% max  x'*A*x - rho*sum(g_p(x_i))
% s.t. x'*Diag(b)*x = 1

% Input:
% A: n*n symmetric PSD matrix
% b: n*1 positive vector, b > 0
% x0: initial point
% rho: regularization parameter, rho >= 0
% type: the surrogate function type, including: 'Lp','log' and 'exp'
% p: 0<p<=1 for 'Lp' or p>0 for 'log' and 'exp'

% Output:
% x: the output
% F_x: the objective values at each iteration

MaxIter = 1000;
tol = 1e-5; % tolerance, stopping criterion

epsi = 1e-8; % smoothing parameter
n = size(A,1);

% Initialize
x = x0;
F_x = zeros(MaxIter,1);

for k = 1:MaxIter
    
    a = A*x;
    
    switch type
        case 'Lp'
            % obj for Lp-norm
            F_x(k) =  x'*a - rho*sum(abs(x).^p);
            % weights for L_p norm
            w = p/2*epsi^(p-2)*ones(n,1);
            ind = (abs(x) > epsi);
            w(ind) = (p/2).*(abs(x(ind)).^(p-2));
        case 'log'
            % obj for log
            F_x(k) = x'*a- rho/log(1+1/p)*sum(log(abs(x)/p+1));
            % weights for log
            w = (1/(2*epsi*(p+epsi)*log(1+1/p)))*ones(n,1);
            ind = (abs(x) > epsi);
            w(ind) = (0.5/log(1+1/p))./(x(ind).^2+p*abs(x(ind)));
        case 'exp'
            % obj for exp
            F_x(k) = x'*a - rho*sum(1-exp(-abs(x)./p));
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
      
    val = min(w./b);
    I_min = (abs(w./b - val)<eps);
    u_min = -rho*val;
    s = sum((b(~I_min).*a(~I_min).^2)./(u_min*b(~I_min) + rho*w(~I_min)).^2);
    
    if any(a(I_min)~=0) || s > 1
        low_bound = max(abs(a)./sqrt(b)-rho*w./b);
        up_bound = sqrt(sum(a.^2./b)) + u_min;
        mu = bisection(a,b,rho*w,low_bound,up_bound);
        x = a./(mu*b+rho*w);
    else
        mu = u_min;
        x(~I_min) = a(~I_min)./(mu*b(~I_min) + rho*w(~I_min));
        x(I_min) = 0;
        ind_last = find(I_min,1,'last');
        x(ind_last) = sqrt((1 - s)/b(ind_last));
    end
    
end


function mu = bisection(a,b,rho_w,low_bound,up_bound)

mid = 0.5*(low_bound+up_bound);
fun_value = sum(b.*(a./(mid*b+rho_w)).^2);
while abs(fun_value-1) > 1e-6
    if fun_value > 1
        low_bound = mid;
    else
        up_bound = mid;
    end
    mid = 0.5*(low_bound+up_bound);
    fun_value = sum(b.*(a./(mid*b+rho_w)).^2);
end
mu = mid;
