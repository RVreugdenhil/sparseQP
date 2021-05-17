function [fobj] = ComputeMainObj(x,A,B,k)
% min_{x} (x'*A*x) / (x'*B*x)
x = proj_l0(x,k);
fobj = (x'*A*x) / (x'*B*x);



