clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal Component Hierarchy for Sparse Quadratic Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simple example of BR and DP on 500 x 500 synthetic dataset

seed = 0;

% Generation of the synthetic data set
N = 500; %amount of rows
n = 500; %amount of columns
sparsity = 10; % nonzero elements

% generate synthetic dataset
    [X,Y,x_true,Itrue] = gendat(N,n,0,0.3,sparsity,6); 
    Q = X'*X/N;%Construct the matrix Q
    c = -2*Y'*X/N;%Construct the vector c

% Implementing our Best Response and Dual Program algorithm for different k

    entire = tic;
    [x_br,~,In_br,~] = BR(Q,c,10, sparsity,20,6,2,4);
    timeElapsed_br = toc(entire);
    % size of the reduced problem
    sizeIn_br = size(In_br,1);
    % in-sample MSE
    [train_br] = evaluateError(x_br, X, Y);


    entire = tic;
    [x_dp,~,In_dp,~] = DP(Q,c,10, sparsity,20,6e-2,12,4,4);
    timeElapsed_dp = toc(entire);
    % size of the reduced problem
    sizeIn_dp = size(In_dp,1);
    % in-sample MSE
    [train_dp] = evaluateError(x_dp, X, Y);
            
