function [d] = lbfgs(g,s,y,Hdiag)
% BFGS Search Direction
%
% This function returns the (L-BFGS) approximate inverse Hessian,
% multiplied by the gradient
% inv(H)*grad
% If you pass in all previous directions/sizes, it will be the same as full BFGS
% If you truncate to the k most recent directions/sizes, it will be L-BFGS
%
% s - previous search directions (n by k)
% y - previous step sizes (n by k)
% g - gradient (n by 1)
% Hdiag - value of initial Hessian diagonal elements (scalar)

[n,d,k] = size(s);

for i = 1:k
    ro(i,1) = 1 / trace(y(:,:,i)'*s(:,:,i));
end

q = zeros(n,d,k+1);
r = zeros(n,d,k+1);
al =zeros(k,1);
be =zeros(k,1);

q(:,:,k+1) = g;

for i = k:-1:1
    al(i) = ro(i) * trace(s(:,:,i)'*q(:,:,i+1));
    q(:,:,i) = q(:,:,i+1) - al(i)*y(:,:,i);
end

% Multiply by Initial Hessian
r(:,:,1) = Hdiag*q(:,:,1);

for i = 1:k
    be(i) = ro(i)*trace(y(:,:,i)'*r(:,:,i));
    r(:,:,i+1) = r(:,:,i) + s(:,:,i)*(al(i)-be(i));
end
d=r(:,:,k+1);

