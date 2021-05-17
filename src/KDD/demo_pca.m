% This program solves the following optimization problem:
% min_x (x'Ax)/(x'Cx), s.t .||x||_0 <= k

clc;clear all;close all;
addpath('util','solver','data');
randn('seed',0); rand('seed',0);



data_id = 6;
k = 15;

[dat_x] = SelectData(data_id);
Cov  = mycov(dat_x);
opt_A  = -Cov;
opt_C = eye(size(Cov,1));
HandleObj = @(x)ComputeMainObj(x,opt_A,opt_C,k);
[x,his,tt] = DEC(opt_A,opt_C,k,[6;6],0);
plot(tt,his)
HandleObj(x)
his(end)



