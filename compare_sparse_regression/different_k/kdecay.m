clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal Component Hierarchy for Sparse Quadratic Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs the sparse regression application for both the BR and DP
% over all possible dimensions of k
%% Data generation of an online source
% set dataset equal to a value between 1 - 5 to select different datasets
dataset = 5;

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

% Get the indices for X in each data set
[x_index,y_index,khat] = get_xy_index(filename);
% Get path to the data folder
fullfilename = strcat(get_path_data(),filename);
A = readmatrix(fullfilename);
Anorm = normalize_column(A(:,[x_index,y_index]));
X = Anorm(:,1:end-1);
Y = Anorm(:,end);

N = height(Y);
n = width(X);
sparsity = 10;
%% Loop over everything
% Different training and testing samples we use
replications =1;
eigtol = 1e-12; %Tolerance of computing the eigenvalues standard 1e-12 
max_iter_br = 40; %T_BR
period_br = 10; %p_BR
period_dp = 100; %p_DP
max_iter_dp = 30; %T_DP
stepsize = 2e-2; % stepsize constant for the DP
ki = [1:n]; % different dimensions of k we test on
eta = sqrt(round(0.7*N)); % eta = sqrt(Ntrain) where 70% of the data is used for testing

% Implementing BR and DP algorithm for different k
for rep = 1:replications   
    seed = rep-1;
    
    % generate training-test data
    [Xtrain,Ytrain,Xtest,Ytest] = traintest(X,Y,seed,0.7);
    Q = (Xtrain'*Xtrain)/size(Xtrain,1);%Construct the matrix Q
    c = (-2*Ytrain'*Xtrain)/size(Xtrain,1);%Construct the vector c^T
    
    % determine big-M constant
    [x_br_lower] = determine_M(Q,c,eta, sparsity,stepsize,20,max_iter_br,period_br);
    M = 4*max(abs(x_br_lower));
    
    % implementing BR and DP for different dimensions of k
    for l = 1:n
        k = ki(l);
        entire = tic;
        [x_br,obj_br(l),In_br,timeloop_br(l,rep)] = BR(Q,c,eta, sparsity,k,max_iter_br,period_br,M);
        % size of the reduced problem
        sizeIn_br(l) = size(In_br,1);
        % in-sample MSE
        [train_br(l,rep)] = evaluateError(x_br, X, Y);
        % [test_br(l)] = evaluateError(x_br, Xtest, Ytest);
        timeElapsed_br(l,rep) = toc(entire);

        entire = tic; 
        [x_dp,obj_dp(l),In_dp,timeloop_dp(l,rep)] = DP(Q,c,eta, sparsity,k,stepsize,max_iter_dp,period_dp,M);
        % size of the reduced problem
        sizeIn_dp(l) = size(In_dp,1);
        % in-sample MSE
        [train_dp(l,rep)] = evaluateError(x_dp, X, Y);
        % [test_dp(l)] = evaluateError(x_dp, Xtest, Ytest);
        timeElapsed_dp(l,rep) = toc(entire);
    end
end
% averaging over all training-test splits
train_br_final = sum(train_br,2)/replications;
train_dp_final = sum(train_dp,2)/replications;
timeloop_br_final = sum(timeElapsed_br,2)/replications;
timeloop_dp_final = sum(timeElapsed_dp,2)/replications;

%%
figure('Renderer', 'painters', 'Position', [10 10 750 300])
hold on
plot(train_br_final,'b', 'LineWidth',2)
plot(train_dp_final,'r', 'LineWidth',2)
lgd = legend('BR','DP')
lgd.FontSize = 14;
xlabel('dimension of subspaces $k$','Interpreter','latex','FontSize', 16) 
ylabel('MSE','FontSize', 16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',15)
axis([0 50 0.012 inf])

figure('Renderer', 'painters', 'Position', [10 10 750 300])
hold on
set(gca, 'YScale', 'log');
plot(timeloop_br_final,'b', 'LineWidth',2)
plot(timeloop_dp_final,'r', 'LineWidth',2)
lgd = legend('BR','DP')
lgd.FontSize = 14;
xlabel('dimension of subspaces $k$','Interpreter','latex','FontSize', 16) 
ylabel('time (seconds)','FontSize', 16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',15)
axis([0 50 10e-2 inf])


