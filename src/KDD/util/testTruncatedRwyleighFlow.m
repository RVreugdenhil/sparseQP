function testTruncatedRwyleighFlow
clc;clear all;close all;
randn('seed',0);
rand('seed',0);
n = 10;
A = randn(n);
A = A*A';
B = randn(n);
B = B*B';

[x] = TruncatedRwyleighFlow(A,B);



function [x] = TruncatedRwyleighFlow(A,B)
% max_x x'Ax/x'Bx, s.t. x(end) = 1

n = size(A,1);
I = eye(n);
x = randn(n,1);
eta = 1/max(abs(eig(B)));

his = [];
for iter =1:1000,
    rho = (x'*A*x)/(x'*B*x);
    his(iter) =rho;
    fprintf('iter:%d, fobj:%f\n',iter,rho);
    C = I + eta/rho * (A-rho*B);
    Cx = C*x;
    x = Cx/norm(Cx);
    x = x/norm(x);
       x(n)=1;
end

plot(his)

