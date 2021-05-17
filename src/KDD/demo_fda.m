% This program solves the following optimization problem:
% min_x (x'Ax)/(x'Cx), s.t .||x||_0 <= k
clc;clear all;close all;
addpath('util','solver','data');
randn('seed',0); rand('seed',0);


global Example
Example = [];
data_id = 6;
k = 8;

[dat_x,dat_y] = SelectData(data_id);
[A,C] = getFisher_AC(dat_x,dat_y);
opt_A = -A;
opt_C = C;
HandleObj = @(x)ComputeMainObj(x,opt_A,opt_C,k);
[x,his,tt] = DEC(opt_A,opt_C,k,[6,6],1);
plot(tt,his)
HandleObj(x)
his(end)
