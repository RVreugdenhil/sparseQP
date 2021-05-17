clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Principal Component Hierarchy for Sparse Quadratic Programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzing the eigenvalue decay and the Frobenius norm ||Q_k - Q||_F
%%
filename = 'Superconductivity.csv';

% Get the indices for X in each data set
[xindex, yindex] = get_xy_index(filename);
% Get path to the data folder
fullfilename = strcat(get_path_data(),filename);
A = readmatrix(fullfilename);
A = normalize_column(A);
X = A(:, xindex);
Y = A(:, yindex);

% Loop over everything
nbReplications = 100;
l4 = zeros(nbReplications,6);
Tol = 1e-12;
k = size(X, 2);

error = zeros(nbReplications, k);
for ita = 1:nbReplications
    seed = ita-1;
    [Xtrain, ~, ~, ~] = traintest(X,Y,seed, 0.7);
    % Compute the matrix Q from Xtrain
    Q = Xtrain'*Xtrain/size(Xtrain, 1);
    
    [V,D] = eigs(Q,k,'largestreal','Tolerance',Tol);    %Calculate the eigenvalues and the eigenvectors of matrix Q
    lambda = diag(D);
    % SelectingLambda
    li = [1,2,3,4,5,10];
    l4(ita,:) = lambda(li);
    
    for k1 = 1:k
        error(ita, k1) = norm(Q - V(:, 1:k1)*diag(lambda(1:k1))*V(:, 1:k1)', 'fro');
    end
end

display(['Ratio lambda_1/lambda_10: ', num2str(mean(l4(:,1)./l4(:,5)))]);


% Plot the eigenvalue decay
for i = 1:length(li)
    [f(i, :),x(i, :)] = ksdensity(l4(:,i)); 
end


figure('Renderer', 'painters', 'Position', [10 10 750 400])
hold on
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
for i = 1:length(li)
    %semilogx(x(i,:), f(i,:));
    loglog(x(i,:), f(i,:), 'LineWidth', 2);
end

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
xlabel('eigenvalue', 'FontSize', 24);
xlim([0.01, 15])
ylabel('density', 'FontSize', 24);
legend('Largest eigevalue','2nd largest eigenvalue','3rd largest eigenvalue','4th largest eigenvalue','5th largest eigenvalue','10th largest eigenvalue', 'FontSize', 14);
% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%%
% Plot reconstruction error for all data sets

allfilename = {'UJIndoor.csv', 'Crime.csv', 'OnlineNewsPopularity.csv', 'Facebook.csv',  'Superconductivity.csv'};
% to compute the percentage of error for hat{k}
percentile = 0.1;
f1 = figure('Renderer', 'painters', 'Position', [10 10 750 400])
hold on

for j = 1:length(allfilename)
    filename = allfilename{j}
    % Get the indices for X in each data set
    [xindex, yindex] = get_xy_index(filename);
    % Get path to the data folder
    fullfilename = strcat(get_path_data(),filename);
    A = readmatrix(fullfilename);
    A = normalize_column(A);
    X = A(:, xindex);
    Y = A(:, yindex);
    %
    % Loop over everything
    nbReplications = 100;
    l4 = zeros(nbReplications,6);
    Tol = 1e-12;
    k = size(X, 2);

    error = zeros(nbReplications, k);
    for ita = 1:nbReplications
        seed = ita-1;
        [Xtrain, ~, ~, ~] = traintest(X,Y,seed, 0.7);
        % Compute the matrix Q from Xtrain
        Q = Xtrain'*Xtrain/size(Xtrain, 1);

        [V,D] = eigs(Q,k,'largestreal','Tolerance',Tol);    %Calculate the eigenvalues and the eigenvalue of matrix Q
        lambda = diag(D);
        % SelectingLambda
        li = [1,2,3,4,5,10];
        l4(ita,:) = lambda(li);

        for k1 = 1:k
            error(ita, k1) = norm(Q - V(:, 1:k1)*diag(lambda(1:k1))*V(:, 1:k1)', 'fro');
        end
        

    end
    indicator = (error < percentile*error(:, 1));
    [~,B]=max(indicator,[],2);
    mean(B)
    
    % plot the mean error
    figure(f1)
    plot(1:k, nanmean(error, 1), 'LineWidth', 2);
end
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
xlabel('Number of subspaces', 'FontSize', 24);
ylabel('Frobenius error', 'FontSize', 24);
legend('UJIndoor $(\hat{k} = 48.87)$', 'Crime $(\hat{k} = 52.13)$', 'OnlineNewsPopularity $(\hat{k} = 20)$', 'Facebook $(\hat{k} = 14)$', 'Superconductivity $(\hat{k} = 9)$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'XScale', 'lin');
set(gca, 'YScale', 'lin');
xlim([0, 100]);
xticks([20 40 60 80 100])
xticklabels({'20','40','60','80','100'})
% The below code expand the axis to minimize the white gap
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];