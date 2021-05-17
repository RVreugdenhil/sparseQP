% This program solves the following optimization problem:
% min_x (x'Ax) / (x'Cx), s.t. ||x||_0<=k

clc;clear all;close all;
addpath('util','solver','data');
randn('seed',1);
rand('seed',1);

[data,ks] = ReadParameters;
HandleFDA = @(opt_A,opt_C,k)DEC(opt_A,opt_C,k,[6;6],1);
times = 5;

for idata = 1:length(data)
    for ik = 1:length(ks)
        cur_data_index = data(idata);
        cur_k = ks(ik);
        for itime = 1:times
            [dat_x,dat_y] = SelectData(cur_data_index);
            [A,B] = getFisher_AC(dat_x,dat_y);
            opt_A = -A;
            opt_C = B;
            [x] = HandleFDA(opt_A,opt_C,cur_k);
            fobj = ComputeMainObj(x,opt_A,opt_C,cur_k);
            fprintf('data: %d, dim: %f, k: %d, fobj: %f\n',cur_data_index,size(A,1),cur_k,fobj);
            objs(itime) = fobj;
        end
        result{idata,ik} = objs;
        result
        save(mfilename,'result');
    end
end









 