function [dataNorm] = normalize_column(A)
% Normalize the matrix A such that each column has range [0, 1]
% Input: 
% A: unnormalized matrix

% Output:
% dataNorm: normalized matrix

    datamin = A - min(A, [], 1); 
    % Scale to the max of each column 
    dataNorm = datamin ./ max(datamin, [], 1); 
end

