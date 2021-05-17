function [x,mu] = acc_steepest_ascent(A,B,x0,T) 

% Compute the leading generalized eigenvector of the matrix pair
% (A,B) using the steepest ascent method with preconditioning and
% accerleration (SQUAREM)

% Input:
% A: n*n symmetric matrix
% B: n*n symmetric positive definite matrix
% x0: initial point
% T: precondition vector, precondtioned with matrix diag(T)

% Output:
% x: the leading generalized eigenvector
% mu: the largest generalized eigenvalue

MaxIter = 1000;
epsilon = sqrt(length(x0))*1e-4; % accurate inner iteration does not accelerate the overall convergence very much
x = x0;
mu = 0;

for k = 1:MaxIter
    % SQUAREM
    [x1,mu,resi_norm] = steepest_ascent_one_step(A,B,x,T);
    
    % Stopping criterion
    if resi_norm <= epsilon
        return;
    end
    
    [x2,mu2,resi_2] = steepest_ascent_one_step(A,B,x1,T);
    r = x1 - x;
    v = x2 - x1 - r;
    alpha = -norm(r)/norm(v);
    
    x_prime = x - 2*alpha*r + alpha^2*v;
    Ax = A*x_prime; Bx = B*x_prime;
    a_xx = x_prime'*Ax; b_xx = x_prime'*Bx;
    mu_prime = a_xx/b_xx;  % eigenvalue

    if mu_prime >= mu % ensure the objective is increasing
        x = x_prime;
    else
        x = x2;
    end   
end

function [x,mu,resi_norm] = steepest_ascent_one_step(A,B,x,T)

Ax = A*x; Bx = B*x;
a_xx = x'*Ax; b_xx = x'*Bx;
mu = a_xx/b_xx;  % eigenvalue
r = Ax - mu*Bx; % residual
resi_norm = norm(r);
r = T.*r; % precondition
Ar = A*r; Br = B*r;
a_rr = r'*Ar; b_rr = r'*Br;
a_rx = x'*Ar; b_rx = x'*Br;

% compute the coefficients of the quadratic polynomial
a = a_rr*b_rx - a_rx*b_rr;
b = a_rr*b_xx - a_xx*b_rr;
c = a_rx*b_xx - a_xx*b_rx;

tao1 = (-b - sqrt(b^2-4*a*c))/(2*a);
tao2 = (-b + sqrt(b^2-4*a*c))/(2*a);
mu1 = (a_xx+2*tao1*a_rx+tao1^2*a_rr)/(b_xx+2*tao1*b_rx+tao1^2*b_rr);
mu2 = (a_xx+2*tao2*a_rx+tao2^2*a_rr)/(b_xx+2*tao2*b_rx+tao2^2*b_rr);
if mu1 > mu2
    tao = tao1;
else
    tao = tao2;
end

% update x
x = x + tao*r;
x = x/sqrt(b_xx+2*tao*b_rx+tao^2*b_rr); % normalize

