clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal Component Hierarchy for Sparse Quadratic Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script runs the sparse regression application for synthetic data
% with different dimensions mentioned in the supplementary of the main paper
%% Loop over everything
measurements =6; %Amount of different values for N we want to test on 
features = 1; %Amount of different values for n we want to test on 

rho = 0.5; %Correlation coefficient
sparsity = 10; %amount of non zero components
SNR = 6; %signal to noise
max_iter_br = 20; %T_BR
period_br = 4; %p_BR
period_ws = 4; % p of the warm start
period_dp = 10; %p_DP
max_iter_dp = 30; % T_DP
stepsize = 3e-2; % stepsize constant DP
max_iter_ws = 20; %Iterations of the warm start
max_iter_Screen = 60; %Iterations of the screening
eta = 10; % regularisation term
M = 4; % compute big M constant
replications = 50; %amount of replications
for tot = 1:replications
seed = tot-1; %seed for reproducability

    for j = 1:measurements
        for i = 1:features

        % Generation of the synthetic data set
        Nj = [100,500,1000,5000,10000,20000];
        N = Nj(j); %amount of rows
        ni = [1000]; 
        n = ni(i); %amount of columns

        [X,Y,~,Itrue] = gendat(N,n,seed,rho,sparsity,SNR); 
        Q = X'*X/height(Y);%Construct the matrix Q
        c = -2*Y'*X/height(Y);%Construct the vector c



        % Implementing BR and DP 
        ki = [30;40;min(N,n)/2];
        for l = 1:length(ki)
            k = ki(l);
            entire = tic;
            [x_br,~,In_br_u,~] = BR(Q,c,eta, sparsity,k,max_iter_br,period_br,M,1e-12);
            tbr(l) = toc(entire);
            % in-sample MSE
            [train_br(l)] = evaluateError(x_br, X, Y);

            
            entire = tic;
            [x_dp,~,In_dp,~] = DP(Q,c,eta, sparsity,k,stepsize,max_iter_dp,period_dp,M,1e-12);
            tdp(l) = toc(entire);
            % in-sample MSE
            [train_dp(l)] = evaluateError(x_dp, X, Y);
 
        end
        % Implementation of the method by Bertsimas & van Parys.
        entire = tic;
        x = zeros(n,1);
        [x_ws,~,~,timeloop_ws] = WarmStart(X,Y,Q,c,eta,sparsity,max_iter_ws,period_ws,M);
        tws = toc(entire);
        % in-sample MSE
        [train_ws] = evaluateError(x_ws, X, Y);
        
        % Implementation of the method by Beck and Eldar
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
        [A,C] = getFisher_AC(X,Y,mean(Y));
        opt_A = -A;
        opt_C = C;
        HandleObj = @(x)ComputeMainObj(x,opt_A,opt_C,sparsity);
        [x,his,tt] = DEC(opt_A,opt_C,sparsity,[6,6],1);
        Ikdd = find(x);
        x = zeros(n,1);
        x(Ikdd) = -0.5*inv(eye(length(Ikdd))/eta+Q(Ikdd,Ikdd))*c(1,Ikdd)';
        tkdd = toc;
        % in-sample MSE
        [train_kdd] = evaluateError(x, X, Y);
    
        % making of the table
        result(20*(i-1)+1:20*i,j,tot) = [train_br';train_dp';train_ws;train_bez;train_bef;train_kdd;...
            tbr';tdp';tws;tbez;tbef;tkdd];
        end    
    end
    tot
end
% Taking the average over all results
resultfinal = sum(result,3)/replications
% 
figure('Renderer', 'painters', 'Position', [10 10 750 300])
hold on
set(gca, 'YScale', 'log');
semilogy(resultfinal([13],:)', 'r', 'LineWidth',2)
semilogy(resultfinal([16],:)', 'b--', 'LineWidth',2)
semilogy(resultfinal([17],:)', 'k-.', 'LineWidth',2)
lgd = legend('DP $k =\min(N,n)/2$','BR $k =\min(N,n)/2$','warm start','Interpreter','latex','Location','northwest')
lgd.FontSize = 14;
xlabel('sample size $N$','Interpreter','latex','FontSize', 16) 
ylabel('time (seconds)','FontSize', 16)
xticks([1 2 3 4 5 6])
xticklabels({'100','500','1000','5000','10000','20000'})
axis([1 inf 10e-2 10e1])

