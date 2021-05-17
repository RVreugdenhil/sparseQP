function [x] = reoptimizer(x,A,B)
% min_{x} x'Ax / x'Bx

[ind]=find(x~=0);
% [V,D]=eig(A(ind,ind),B(ind,ind));
% x1 = V(:,1);
sigma = 'sa';
[x1,~]=eigs(A(ind,ind),B(ind,ind),1,sigma);
x1(x1==0)=eps;
x(ind) = x1;

