clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal Component Hierarchy for Sparse Quadratic Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs the sparse regression application for synthetic data
% of constant dimension with different variables rho,eta,s and SNR

replications = 50;
variables = 5;
% Generation of the synthetic data set
N = 1000; %amount of rows
n = 1000; %amount of columns
max_iter_br = 20; %T_BR
max_iter_dp = 500; %T_DP
stepsize = 4e-3; % stepsize constant of the DP method
max_iter_ws = 20; %Iterations of the method by Bertsimas & van Parys
period_br = 4; % p_BR
period_ws = 4; % p of the warm start method
period_dp = 10; % p_DP
max_iter_Screen = 60; %Iterations of the method by Gomez & Atamturk
eigtol = 1e-12; %Tolerance of computing the eigenvalues standard 1e-12 
ki = [20;30;40;200;400;500;600;800]; 
etai = [100;10;1;0.1;0.01;0.001];
si = [5;10;20;30;40];% differnt values of s we test on
rhoi = [0.7;0.8;0.9]; % different values rho we test on
SNRi = [20;6;3;1;0.05];% different SNR we test on
M = 4; % determine the big-M constant
% define which method we want to test on 
% met = 1 --> SNR
% met = 2 --> RHO
% met = 3 --> S
% met = 4 --> eta
%for met = 1:4
met = 2;
for rep = 1:replications
    seed = rep-1; %seed for reproducability
    
    if met ==4
        variables = 6;
    elseif met == 2
        variables = 3;
    end
    
    
    for i = 1:variables
        
        % determining which variable we test on
        if met ==1
            SNR = SNRi(i);
        else
            SNR = 6;
        end
        
        if met ==2
            rho = rhoi(i);
        else
            rho = 0.5;
        end

        if met ==3
            sparsity = si(i);
        else
            sparsity = 10;
        end
        
        if met==4
            eta = etai(i);
        else
            eta = 10;
        end
        
        % generate synthetic dataset
        [X,Y,~,~] = gendat(N,n,seed,rho,sparsity,SNR); 
        Q = X'*X/N;%Construct the matrix Q
        c = -2*Y'*X/N;%Construct the vector c

       % Implementing our Best Response and Dual Program algorithm for different k
        
        for l = 1:length(ki)
            k = ki(l);
            entire = tic;
            [x_br,~,In_br,~] = BR(Q,c,eta, sparsity,k,max_iter_br,period_br,M);
            tbr(l) = toc(entire);
            % size of the reduced problem
            sizeIn_br(l) = size(In_br,1);
            % in-sample MSE
            [train_br(l)] = evaluateError(x_br, X, Y);
            
          
            entire = tic;
            [x_dp,~,In_dp,~] = DP(Q,c,eta, sparsity,k,stepsize,max_iter_dp,period_dp,M);
            tdp(l) = toc(entire);
            % size of the reduced problem
            sizeIn_dp(l) = size(In_dp,1);
            % in-sample MSE
            [train_dp(l)] = evaluateError(x_dp, X, Y);
            
        end
        
        % Implementation of the method by Bertsimas & van Parys.
        entire = tic;
        [x_ws,~,~,~] = WarmStart(X,Y,Q,c,eta,sparsity,max_iter_ws,period_br,M);
        tws = toc(entire);
        % in-sample MSE
        [train_ws] = evaluateError(x_ws, X, Y);
        
        
        % Implementation of the method by Gomez & Atamturk
        
        % only test for different eta
        if met==4
        entire = tic;
        [x_screen,~,~,sizeIn_screen] = Screening(X,Y,Q,c,eta,sparsity,max_iter_Screen,M);
        tscreen = toc(entire);
        % in-sample MSE
        [train_screen] = evaluateError(x_screen, X, Y);
        else
        sizeIn_screen = [];
        train_screen= [];
        timeloop_screen = [];
        tscreen = [];
        end
        
        % Implementation of the method by Beck and Eldar
        x = zeros(n,1);
        tic
        [x_bez, ~] = beckzero(Q,c,x,eta,sparsity);
        tbez = toc;
        % in-sample MSE
        [train_bez] = evaluateError(x_bez, X, Y);

        tic
        [x_bef, ~, I_f] = beckfull(Q,c,x,eta,sparsity);
        tbef = toc;
        % in-sample MSE
        [train_bef] = evaluateError(x_bef, X, Y);
        
        if train_bef > train_bez
            train_bef = train_bez;
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
        [train_kdd] = evaluateError(x_kdd, X, Y);
    

        % summarize all data in one table
        result(:,i,rep) = [sizeIn_dp';sizeIn_br';sizeIn_screen;...
            train_dp';train_br';train_ws;train_screen;train_bez;train_bef;train_kdd;...
            tdp';tbr';tws;tscreen;tbez;tbef;tkdd];

    end
rep
end
% Taking the average over all results
resultfinal = sum(result,3)/replications
