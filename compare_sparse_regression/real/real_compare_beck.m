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
for dat = 1:5
dataset = dat;

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
replications =50;
eigtol = 1e-12; %Tolerance of computing the eigenvalues standard 1e-12 
max_iter_br = 40; %T_BR
period_br = 10; %p_BR
period_dp = 100; %p_DP
max_iter_dp = 3000; %T_DP
stepsize = 2e-3; %stepsize constant in the dual program
max_iter_ws = 20; %Iterations of the method by Bertsimas & van Parys
max_iter_Screen = 60; %Iterations of the method by Gomez & Atamturk

%different values for k we test on
ki = [khat];
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
%             entire = tic;
%             [x_br,obj_br(l),In_br_u,~] = BR(Q,c,eta, sparsity,k, max_iter_br,period_br,M);
%             % size of the reduced problem
%             sizeIn_br(l) = size(In_br_u,1);
%             % in-sample MSE
%             [train_br(l)] = evaluateError(x_br, Xtrain, Ytrain);
%             % out-sample MSE
%             [test_br(l)] = evaluateError(x_br, Xtest, Ytest);
%             timeElapsed_br(l) = toc(entire);
          
            entire = tic;
            [x_dp,obj_dp(l),In_dp,~] = DP(Q,c,eta, sparsity,k,stepsize,max_iter_dp,period_dp,M);
            % size of the reduced problem
            sizeIn_dp(l) = size(In_dp,1);
            % in-sample MSE
            [train_dp(l)] = evaluateError(x_dp, Xtrain, Ytrain);
            % out-sample MSE
            [test_dp(l)] = evaluateError(x_dp, Xtest, Ytest);
            timeElapsed_dp(l) = toc(entire);
        end
        avgY = sum(Ytrain)/length(Ytrain);
        tic
        [A,C] = getFisher_AC(Xtrain,Ytrain,avgY);
        opt_A = -A;
        opt_C = C;
        HandleObj = @(x)ComputeMainObj(x,opt_A,opt_C,sparsity);
        [x,his,tt] = DEC(opt_A,opt_C,sparsity,[6,6],1);
        Ikdd = find(x);
        x = zeros(n,1);
        x(Ikdd) = -0.5*inv(eye(length(Ikdd))/eta+Q(Ikdd,Ikdd))*c(1,Ikdd)';
        tkdd = toc;
        %         plot(tt,his)
%         HandleObj(x)
%         his(end)
        [train_kdd] = evaluateError(x, Xtrain, Ytrain);
        [test_kdd] = evaluateError(x, Xtest, Ytest);
%         x = zeros(n,1);
%         tic
%         [x, ~] = beckzero(Q,c,x,eta,sparsity);
%         tzero = toc;
%         [train_z] = evaluateError(x, Xtrain, Ytrain);
%         [test_z] = evaluateError(x, Xtest, Ytest);
%         %Iz = find(x);
%         x = zeros(n,1);
%         tic
%         [x, ~, I_f] = beckfull(Q,c,x,eta,sparsity);
%         tfull = toc;
%         [train_f] = evaluateError(x, Xtrain, Ytrain);
%         [test_f] = evaluateError(x, Xtest, Ytest);
%         
%         if train_f > train_z
%             train_f = train_z;
%             test_f = test_z;
%         end
        
        % Implementation of the method by Bertsimas & van Parys.
%         if dataset == 4 %else this method will run out of memory
%             size_ws = [];
%             obj_ws = [];
%             timeElapsed_ws = [];
%             train_ws = [];
%             test_ws = [];
%         else
%         entire = tic;
%         [x_ws,obj_ws,~,~] = WarmStart(Xtrain,Ytrain,Q,c,eta,sparsity,max_iter_ws,period_br,M);
%         % in-sample MSE
%         [train_ws] = evaluateError(x_ws, Xtrain, Ytrain);
%         % out-sample MSE
%         [test_ws] = evaluateError(x_ws, Xtest, Ytest);
%         timeElapsed_ws = toc(entire);
%         end
%         % Implementation of the method by Gomez & Atamturk
%         entire = tic;
%         [x_screen,obj_screen,~,sizeIn_screen] = Screening(Xtrain,Ytrain,Q,c,eta,sparsity,max_iter_Screen,M);
%         % in-sample MSE
%         [train_screen] = evaluateError(x_screen, Xtrain, Ytrain);
%         % out-sample MSE
%         [test_screen] = evaluateError(x_screen, Xtest, Ytest);
%         timeElapsed_screen = toc(entire);

        % makeing a table with all the results
            result(:,j) = [train_dp';train_kdd;...
                test_dp';test_kdd;...
                timeElapsed_dp';tkdd];%...
j
end
% averaging over all different training-test split
resultfinal = sum(result,2)/replications

%% saving the different results
if dataset == 1 
%     save('OnlineNewsPopularity.mat','result')
    xlswrite('OnlineNewsPopularity_kdd',resultfinal);
% %     standard = std(resultfinal,[],3,'omitnan');
%      xlswrite('OnlineNewsPopularity_beck',standard);
elseif dataset == 2
%     save('Crime.mat','result')
    xlswrite('Crime_kdd',resultfinal);
%     standard = std(resultfinal,[],3,'omitnan');
%     xlswrite('Crime_beck',standard);
elseif dataset == 3
%     save('UJIndoor.mat','result')
    xlswrite('UJIndoor_kdd',resultfinal);
%     standard = std(resultfinal,[],3,'omitnan');
%     xlswrite('UJIndoor_beck',standard);
elseif dataset == 4
%     save('Facebook.mat','result')
    xlswrite('Facebook_kdd',resultfinal);
%     standard = std(resultfinal,[],3,'omitnan');
%     xlswrite('Facebook_beck',standard);
elseif dataset == 5
%     save('Superconductivity.mat','result')
    xlswrite('Superconductivity_kdd',resultfinal);
%     standard = std(resultfinal,[],3,'omitnan');
%     xlswrite('Superconductivity_beck',standard);
end
end