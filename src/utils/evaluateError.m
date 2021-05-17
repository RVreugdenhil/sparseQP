function [error] = evaluateError(w, X, Y)
% Return the average residual
% Input:
% w: regressor (n x 1)
% X: matrix (N x n)
% Y: vector (N x 1)
    N = length(Y);
    error = norm(X*w - Y, 2)^2/N;
end

