function [Xtrain,Ytrain,Xtest,Ytest] = traintest(X,Y,seed, split)
% Create a random split of the data (X,Y) into training and test set
% Input:
% (X, Y): input data
% seed: random seed number
% split: percentage of training between (0, 1)

% Output:
% Xtrain, Ytrain: training set 
% Xtest, Ytest: test set

    N = height(Y);
    Nnew = round(split*N); %Select training data
    rng(seed,'v5uniform') 
    
    Xtrue = randperm(N,Nnew); %Select random indices that form our training set
    Xtest = X;
    Xtest(Xtrue,:) = []; %Delete the training set out of the testing set
    Ytest = Y;
    Ytest(Xtrue) = []; %Delete the training set out of the testing set
    Xtrain = X(Xtrue,:); %Only select the training set
    Ytrain = Y(Xtrue); %Only select the training set
end