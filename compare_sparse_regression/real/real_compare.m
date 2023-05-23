clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal Component Hierarchy for Sparse Quadratic Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs the sparse regression application for real data sets
% mentioned in the paper

%% Data generation of an online source
% set dataset equal to a value between 1 - 5 to select different datasets
dataset = 2;

if dataset == 1 
    filename = 'OnlineNewsPopularity.csv';
elseif dataset == 2
    filename = 'Crime.csv';
elseif dataset == 3
    filename = 'UJIndoor.csv';
elseif dataset == 4
    filename = 'Facebook.csv';
elseif dataset == 5
    filename = 'Superconductivity.csv';
end
result = [];

% Get the indices for X in each data set
[x_index,y_index,khat] = get_xy_index(filename);
% Get path to the data folder
fullfilename = strcat(get_path_data(),filename);
A = readmatrix(fullfilename);
Anorm = normalize_column(A(:,[x_index,y_index]));
X = Anorm(:,1:end-1);
Y = Anorm(:,end);

% Dimensions of the matrix
N = size(X,1);
n = size(X,2);
sparsity = 10;
%% Loop over everything

% Different training and testing splits we use
replications =25;
eigtol = 1e-12; %Tolerance of computing the eigenvalues standard 1e-12 
max_iter_br = 40; %T_BR
period_br = 10; %p_BR
period_dp = 100; %p_DP
max_iter_dp = 3000; %T_DP
stepsize = 2e-3; %stepsize constant in the dual program
max_iter_ws = 20; %Iterations of the method by Bertsimas & van Parys
max_iter_Screen = 60; %Iterations of the method by Gomez & Atamturk

%different values for k we test on
ki = [20,30,40,khat];
split = 0.7;
eta = sqrt(round(split*N));
for j = 1:replications
    
    % Generation of the training and testing data set
    seed = j-1;
    [Xtrain,Ytrain,Xtest,Ytest] = traintest(X,Y,seed,split);
    
    %Construct the matrix Q
    Q = (Xtrain'*Xtrain)/size(Xtrain,1);
    %Construct the vector c^T
    c = (-2*Ytrain'*Xtrain)/size(Xtrain,1);
    
    % compute the big M constant
    [x_inf] = determine_M(Q,c,eta, sparsity,stepsize, 40, max_iter_dp,period_dp);
    M = 4*max(abs(x_inf));
    
    
        % Implementing our best response and dual program algorithm for different k
        
        for l = 1:length(ki)
            k = ki(l);
            entire = tic;
            [x_br,obj_br(l),In_br_u,~] = BR(Q,c,eta, sparsity,k, max_iter_br,period_br,M);
            % size of the reduced problem
            sizeIn_br(l) = size(In_br_u,1);
            % in-sample MSE
            [train_br(l)] = evaluateError(x_br, Xtrain, Ytrain);
            % out-sample MSE
            [test_br(l)] = evaluateError(x_br, Xtest, Ytest);
            tbr(l) = toc(entire);
          
            entire = tic;
            [x_dp,obj_dp(l),In_dp,~] = DP(Q,c,eta, sparsity,k,stepsize,max_iter_dp,period_dp,M);
            % size of the reduced problem
            sizeIn_dp(l) = size(In_dp,1);
            % in-sample MSE
            [train_dp(l)] = evaluateError(x_dp, Xtrain, Ytrain);
            % out-sample MSE
            [test_dp(l)] = evaluateError(x_dp, Xtest, Ytest);
            tdp(l) = toc(entire);
        end
        % Implementation of the method by Bertsimas & van Parys.
        if dataset == 4 %else this method will run out of memory
            size_ws = [];
            obj_ws = [];
            tws = [];
            train_ws = [];
            test_ws = [];
        else
%         entire = tic;
%         [x_ws,obj_ws,~,~] = WarmStart(Xtrain,Ytrain,Q,c,eta,sparsity,max_iter_ws,period_br,M);
%         % in-sample MSE
%         [train_ws] = evaluateError(x_ws, Xtrain, Ytrain);
%         % out-sample MSE
%         [test_ws] = evaluateError(x_ws, Xtest, Ytest);
%         tws = toc(entire);
%         end
%         % Implementation of the method by Gomez & Atamturk
%         entire = tic;
%         [x_screen,obj_screen,~,sizeIn_screen] = Screening(Xtrain,Ytrain,Q,c,eta,sparsity,max_iter_Screen,M);
%         % in-sample MSE
%         [train_screen] = evaluateError(x_screen, Xtrain, Ytrain);
%         % out-sample MSE
%         [test_screen] = evaluateError(x_screen, Xtest, Ytest);
%         tscreen = toc(entire);
        
        % Implementation of the method by Beck and Eldar
        x = zeros(n,1);
        tic
        [x_bez, ~] = beckzero(Q,c,x,eta,sparsity);
        tbez = toc;
        % in-sample MSE
        [train_bez] = evaluateError(x_bez, Xtrain, Ytrain);
        % out-sample MSE
        [test_bez] = evaluateError(x_bez, Xtest, Ytest);

        tic
        [x_bef, ~, I_f] = beckfull(Q,c,x,eta,sparsity);
        tbef = toc;
        % in-sample MSE
        [train_bef] = evaluateError(x_bef, Xtrain, Ytrain);
        % out-sample MSE
        [test_bef] = evaluateError(x_bef, Xtest, Ytest);
        
        if train_bef > train_bez
            train_bef = train_bez;
            test_bef = test_bez;
        end

        % Method of Yuan et al. copied from the authors Github
        tic
        Alg = @(sparsity,c,Q,x0)BlockDecAlg_c(eps*randn(size(x0)),Q,c,sparsity,[12;6],0);
        x0 = randn(size(Q,2),1); 
        HandleObj = @(x)0.5*norm(Q*proj_l0(x,k)-c,'fro')^2;
        [x,his]  = Alg8(sparsity,c,Q,x0); 
        Ikdd = find(x);
        x_kdd = zeros(n,1);
        x_kdd(Ikdd) = -0.5*inv(eye(length(Ikdd))/eta+Q(Ikdd,Ikdd))*c(1,Ikdd)';
        tkdd = toc;
        % in-sample MSE
        [train_kdd] = evaluateError(x_kdd, Xtrain, Ytrain);
        % out-sample MSE
        [test_kdd] = evaluateError(x_kdd, Xtest, Ytest);

        % makeing a table with all the results
            result(:,j) = [sizeIn_dp';sizeIn_br';sizeIn_screen;...
                train_dp';train_br';train_ws;train_screen;train_bez;train_bef;train_kdd;...
                test_dp';test_br';test_ws;test_screen;test_bez;test_bef;test_kdd;...
                tdp';tbr';tws;tscreen;tbez;tbef;tkdd];
j
end
% averaging over all different training-test split
resultfinal = sum(result,2)/replications
