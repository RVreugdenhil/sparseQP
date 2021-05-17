% This program solves the following optimization problem:
% min_x (x'Ax) / (x'Cx), s.t. ||x||_0<=k

clc;clear all;close all;
addpath('util','solver','data');
randn('seed',0);
rand('seed',0);

[data,ks] = ReadParameters;
HandleFDA = @(A,B,k)SGkPCA(A,B,k,3);
times = 1;

for idata = 1:length(data)
    for ik = 1:length(ks)
        cur_data_index = data(idata);
        cur_k = ks(ik);
        for itime = 1:times
            [dat_x,dat_y] = SelectData(cur_data_index);
            [A,B] = getFisher_AC(dat_x,dat_y);
            size(A)
            size(dat_x)
            continue;
            opt_A = -A;
            opt_C = B;
            [x] = HandleFDA(A,B,cur_k);
            fobj = ComputeMainObj(x,opt_A,opt_C,cur_k);
            fprintf('data: %d, dim: %f, k: %d, fobj: %f\n',cur_data_index,size(A,1),cur_k,fobj);
            objs(itime) = fobj;
        end
%         result{idata,ik} = objs;
%         result
%         save(mfilename,'result');
    end
end








 

