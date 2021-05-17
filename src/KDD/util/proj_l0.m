function [x] = proj_l0(x,k)
% min_{X}  0.5 ||X-A||_F^2, s.t. ||X||_0<=k
% x can be a matrix
absx = abs(x);
[~,val] = sort(absx(:),'descend');
zero_ind = val((k+1):end);
x(zero_ind)=0;

